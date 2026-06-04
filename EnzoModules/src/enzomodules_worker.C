// ADR-0005 #3 — the C++ bridge worker process.
//
// A standalone (non-Julia) process that hosts the live Enzo hierarchy and serves
// the EnzoNG bridge over the SAME wire protocol as the Julia reference worker
// (lib/EnzoLib/src/rpc.jl): a line-based control channel on stdin/stdout, and a
// shared file for bulk array arguments.  Because it carries no Julia runtime, its
// gcc/libstdc++ C++ stack never meets Julia's libc++ — the collision that blocked
// in-process MPI (ADR-0004 / ADR-0005) simply cannot occur here.  This is the
// serial worker; the MPI build adds MPI_Init around the same loop (#3b) and
// mpiexec launches N ranks (#4).
//
// The per-symbol typed dispatch is GENERATED from the bridge manifest
// (tools/gen_worker_dispatch.jl → enzomodules_worker_dispatch.inc); the contract
// hash baked there is presented at the handshake so a stale-regeneration mismatch
// is refused, not run.
//
//   usage:  enzomodules_worker <shm_path> <bridge_dylib_path>
//
// Wire protocol (must match rpc.jl exactly):
//   handshake : "READY <hash-decimal>\n"
//   request   : "CALL <symbol> <tok>...\n"   |   "QUIT\n"
//   token     : i<dec> | f<hexbits> | p<dec> | s<base64> | b<off>,<len>,<tag>
//   reply     : "RET <tok>\n"  (tok: i.. / f.. / p.. / "void")   |   "ERR <msg>\n"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <dlfcn.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

// ── a decoded argument (scalar value, or a pointer into the mmap'd shm) ───────
struct Arg {
    char kind = 0;        // 'i' 'f' 'p' 's' 'b'
    long long i = 0;      // scalar int
    double d = 0.0;       // scalar double
    void* p = nullptr;    // ptr scalar OR buffer base+offset (set after mmap)
    std::string s;        // decoded cstring
    size_t off = 0, len = 0;  // buffer descriptor
    char btag = 0;            // buffer eltype tag: d/i/l
};

// ── reply encoders (token only; the loop prints the "RET " prefix) ────────────
static void reply_void(std::string& out)              { out = "void"; }
static void reply_int(std::string& out, int r)        { char b[32]; snprintf(b, sizeof b, "i%d", r); out = b; }
static void reply_ptr(std::string& out, void* r)      { char b[32]; snprintf(b, sizeof b, "p%llu", (unsigned long long)(uintptr_t)r); out = b; }
static void reply_double(std::string& out, double r)  {
    uint64_t u; std::memcpy(&u, &r, sizeof u);
    char b[32]; snprintf(b, sizeof b, "f%llx", (unsigned long long)u); out = b;
}

#include "enzomodules_worker_dispatch.inc"   // GENERATED: worker_dispatch(...) + WORKER_CONTRACT_HASH

// ── base64 decode (only used for Cstring args, e.g. the problem-file path) ────
static std::string b64decode(const std::string& in) {
    static const std::string T =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    int rev[256]; for (int i = 0; i < 256; ++i) rev[i] = -1;
    for (int i = 0; i < 64; ++i) rev[(unsigned char)T[i]] = i;
    std::string out; int val = 0, bits = -8;
    for (unsigned char c : in) {
        if (c == '=' || rev[c] == -1) break;
        val = (val << 6) | rev[c]; bits += 6;
        if (bits >= 0) { out.push_back(char((val >> bits) & 0xFF)); bits -= 8; }
    }
    return out;
}

// ── shm: open + mmap the whole shared file (the worker maps it per-call so it
// always sees the size the client just grew it to; MAP_SHARED makes the kernel's
// OUT-buffer writes visible to the client's subsequent seek/read). ────────────
struct Shm { int fd = -1; void* base = MAP_FAILED; size_t size = 0; };
static bool shm_map(const char* path, Shm& m) {
    m.fd = open(path, O_RDWR);
    if (m.fd < 0) return false;
    struct stat st; if (fstat(m.fd, &st) != 0) { close(m.fd); return false; }
    m.size = (size_t)st.st_size;
    if (m.size == 0) { m.base = nullptr; return true; }   // no buffers this call
    m.base = mmap(nullptr, m.size, PROT_READ | PROT_WRITE, MAP_SHARED, m.fd, 0);
    if (m.base == MAP_FAILED) { close(m.fd); return false; }
    return true;
}
static void shm_unmap(Shm& m) {
    if (m.base && m.base != MAP_FAILED && m.size) { msync(m.base, m.size, MS_SYNC); munmap(m.base, m.size); }
    if (m.fd >= 0) close(m.fd);
}

// ── parse one whitespace-free token into an Arg ───────────────────────────────
static bool parse_token(const std::string& t, Arg& a) {
    if (t.empty()) return false;
    a.kind = t[0];
    const char* rest = t.c_str() + 1;
    switch (a.kind) {
        case 'i': a.i = strtoll(rest, nullptr, 10); return true;
        case 'p': a.p = (void*)(uintptr_t)strtoull(rest, nullptr, 10); return true;
        case 'f': { uint64_t u = strtoull(rest, nullptr, 16); std::memcpy(&a.d, &u, sizeof u); return true; }
        case 's': a.s = b64decode(rest); return true;
        case 'b': {  // b<off>,<len>,<tag>
            char tag = 0; long long off = 0, len = 0;
            if (sscanf(rest, "%lld,%lld,%c", &off, &len, &tag) != 3) return false;
            a.off = (size_t)off; a.len = (size_t)len; a.btag = tag; return true;
        }
        default: return false;
    }
}

#ifdef USE_MPI
// Each rank's OWN session handle.  session_init runs on every rank and returns a
// DIFFERENT pointer per rank (each manages that rank's local grids), but the client
// only ever holds rank 0's.  So when a broadcast call carries the handle, every rank
// substitutes its own; and we capture it from any handle-returning call.  The handle
// is the only scalar pointer ('p') arg in the bridge — all other pointers are shm
// buffers ('b') — so a 'p' token is unambiguously the session handle.
static void* g_my_handle = nullptr;
#endif

// Execute one control line collectively (every rank runs this for the SAME line,
// so collective bridge calls — session_init→CommunicationPartitionGrid, set_boundary,
// compute_dt's CommunicationMinValue, update_from_finer — stay in lockstep).  Writes
// the full reply line ("RET …"/"ERR …") into `reply`.  Returns false on QUIT.
static bool process_command(void* h, const char* shm_path, const std::string& line, std::string& reply) {
    std::istringstream ss(line);
    std::string cmd; ss >> cmd;
    if (cmd == "QUIT") return false;
    if (cmd != "CALL") { reply = "ERR bad-command " + cmd; return true; }

    std::string sym; ss >> sym;
    std::vector<Arg> args;
    std::string tok; bool ok = true; bool has_buf = false;
    while (ss >> tok) {
        Arg a;
        if (!parse_token(tok, a)) { ok = false; break; }
        if (a.kind == 'b') has_buf = true;
        args.push_back(std::move(a));
    }
    if (!ok) { reply = "ERR bad-token in " + sym; return true; }

    // map shm and point each buffer arg into it (IN bytes already present).  NOTE
    // (multi-rank): buffers reference rank-0's shared file; cross-rank field RPC is a
    // #4 concern (collective evolve + a global reduction return only scalars).
    Shm m;
    if (has_buf) {
        if (!shm_map(shm_path, m)) { reply = "ERR shm-map-failed " + sym; return true; }
        for (auto& a : args)
            if (a.kind == 'b') a.p = (void*)((char*)m.base + a.off);
    }
#ifdef USE_MPI
    for (auto& a : args) if (a.kind == 'p') a.p = g_my_handle;   // use THIS rank's handle
#endif

    std::string out;
    void* fn = dlsym(h, sym.c_str());
    bool handled = fn ? worker_dispatch(sym, fn, args, out) : false;
    if (has_buf) shm_unmap(m);   // msync flushes OUT-buffer writes back to the file

    reply = !fn ? ("ERR unknown-symbol " + sym)
          : !handled ? ("ERR undispatched " + sym)
          : ("RET " + out);
#ifdef USE_MPI
    if (reply.rfind("RET p", 0) == 0)                            // captured a new handle
        g_my_handle = (void*)(uintptr_t)strtoull(reply.c_str() + 5, nullptr, 10);
#endif
    return true;
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s <shm_path> <bridge_dylib>\n", argv[0]); return 2; }
    const char* shm_path = argv[1];
    const char* lib_path = argv[2];

    int rank = 0;
#ifdef USE_MPI
    // The worker OWNS MPI in its own process (no Julia runtime here → no C++ ABI
    // collision).  CommunicationInitialize in the bridge is guarded on
    // MPI_Initialized, so it picks up this world instead of re-initializing.
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    // CLAIM stdout for the protocol: the hosted Enzo library prints diagnostics to
    // stdout (e.g. "MPI_Init: NumberOfProcessors = 2"), which would corrupt the
    // control channel — and under mpiexec EVERY rank's stdout is merged into the
    // one the client reads.  So on every rank redirect fd 1 → stderr (noise stays
    // visible for debugging, off the channel); rank 0 keeps a private dup of the
    // real stdout as the control FILE.  Do this BEFORE dlopen (static-init prints).
    FILE* ctrl = nullptr;
    if (rank == 0) ctrl = fdopen(dup(STDOUT_FILENO), "w");
    dup2(STDERR_FILENO, STDOUT_FILENO);
    auto emit = [&](const std::string& s) { if (ctrl) { fputs(s.c_str(), ctrl); fputc('\n', ctrl); fflush(ctrl); } };

    void* h = dlopen(lib_path, RTLD_NOW | RTLD_GLOBAL);
    if (!h) {
        fprintf(stderr, "worker[rank %d]: dlopen(%s) failed: %s\n", rank, lib_path, dlerror());
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 3);
#endif
        return 3;
    }

    // Only rank 0 speaks the control channel; non-zero ranks execute broadcast
    // commands and stay silent on the channel.
    {
        char b[64]; snprintf(b, sizeof b, "READY %llu", (unsigned long long)WORKER_CONTRACT_HASH);
        emit(b);
    }

#ifndef USE_MPI
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        std::string reply;
        if (!process_command(h, shm_path, line, reply)) break;
        emit(reply);
    }
#else
    // Master-driven SPMD: rank 0 reads a command and broadcasts it to all ranks;
    // every rank executes it (collective bridge calls run in lockstep); rank 0 alone
    // replies.  Empty lines are skipped on rank 0 before the broadcast.
    for (;;) {
        std::string line;
        int len = -1;   // -1 = EOF/QUIT sentinel
        if (rank == 0) {
            while (std::getline(std::cin, line) && line.empty()) { /* skip blanks */ }
            len = (std::cin.good() || !line.empty()) ? (int)line.size() : -1;
        }
        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (len < 0) break;                       // EOF on rank 0
        line.resize(len);
        if (len > 0) MPI_Bcast(&line[0], len, MPI_CHAR, 0, MPI_COMM_WORLD);

        std::string reply;
        bool cont = process_command(h, shm_path, line, reply);
        int go = cont ? 1 : 0;
        MPI_Bcast(&go, 1, MPI_INT, 0, MPI_COMM_WORLD);   // agree on QUIT across ranks
        if (!go) break;
        emit(reply);                              // no-op on non-zero ranks (ctrl==null)
    }
#endif

    dlclose(h);
#ifdef USE_MPI
    MPI_Finalize();
#endif
    if (ctrl) fclose(ctrl);
    return 0;
}
