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

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s <shm_path> <bridge_dylib>\n", argv[0]); return 2; }
    const char* shm_path = argv[1];
    const char* lib_path = argv[2];

    void* h = dlopen(lib_path, RTLD_NOW | RTLD_LOCAL);
    if (!h) { fprintf(stderr, "worker: dlopen(%s) failed: %s\n", lib_path, dlerror()); return 3; }

    // handshake: present the baked contract hash in decimal (client parses base 10).
    std::cout << "READY " << (unsigned long long)WORKER_CONTRACT_HASH << "\n" << std::flush;

    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string cmd; ss >> cmd;
        if (cmd == "QUIT") break;
        if (cmd != "CALL") { std::cout << "ERR bad-command " << cmd << "\n" << std::flush; continue; }

        std::string sym; ss >> sym;
        std::vector<Arg> args;
        std::string tok; bool ok = true; bool has_buf = false;
        while (ss >> tok) {
            Arg a;
            if (!parse_token(tok, a)) { ok = false; break; }
            if (a.kind == 'b') has_buf = true;
            args.push_back(std::move(a));
        }
        if (!ok) { std::cout << "ERR bad-token in " << sym << "\n" << std::flush; continue; }

        // map shm and point each buffer arg into it (IN bytes already present).
        Shm m;
        if (has_buf) {
            if (!shm_map(shm_path, m)) { std::cout << "ERR shm-map-failed " << sym << "\n" << std::flush; continue; }
            for (auto& a : args)
                if (a.kind == 'b') a.p = (void*)((char*)m.base + a.off);
        }

        std::string out;
        void* fn = dlsym(h, sym.c_str());
        bool handled = false;
        if (!fn) {
            out.clear();
        } else {
            handled = worker_dispatch(sym, fn, args, out);   // GENERATED: writes OUT bufs into shm
        }
        if (has_buf) shm_unmap(m);   // msync flushes OUT-buffer writes back to the file

        if (!fn)            std::cout << "ERR unknown-symbol " << sym << "\n" << std::flush;
        else if (!handled)  std::cout << "ERR undispatched " << sym << "\n" << std::flush;
        else                std::cout << "RET " << out << "\n" << std::flush;
    }
    dlclose(h);
    return 0;
}
