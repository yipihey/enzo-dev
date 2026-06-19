#!/usr/bin/env julia
# Generate the C++ worker's dispatch table from the bridge manifest (ADR-0005 #3).
#
# The manifest is parsed out of session.jl's @xcall sites (EnzoLib.manifest()), so
# this generator and the Julia client share ONE source of truth for the bridge
# surface.  The emitted file gives the C++ worker, for every C symbol, a typed
# `dlsym`+call with the exact signature — the one thing C++ can't do reflectively.
# The contract hash is BAKED in, so a stale regeneration (worker built from a
# different session.jl than the client) is caught at the handshake, not silently.
#
# Usage:  julia --project=<EnzoLib/test> gen_worker_dispatch.jl [out.inc]
#         (defaults to EnzoModules/src/enzomodules_worker_dispatch.inc)

using EnzoLib

# Julia bridge type  →  (C signature type, C argument expression at index k)
const CTYPE = Dict(
    "Cint"           => "int",
    "Cdouble"        => "double",
    "Cstring"        => "const char*",
    "Handle"         => "void*",
    "Ptr{Cdouble}"   => "double*",
    "Ptr{Cint}"      => "int*",
    "Ptr{Clonglong}" => "long long*",
)
argexpr(ty::String, k::Int) =
    ty == "Cint"           ? "(int)a[$k].i"          :
    ty == "Cdouble"        ? "a[$k].d"               :
    ty == "Cstring"        ? "a[$k].s.c_str()"       :
    ty == "Handle"         ? "a[$k].p"               :
    ty == "Ptr{Cdouble}"   ? "(double*)a[$k].p"      :
    ty == "Ptr{Cint}"      ? "(int*)a[$k].p"         :
    ty == "Ptr{Clonglong}" ? "(long long*)a[$k].p"   :
    error("no C arg-expr for bridge type $ty")

argtypes_of(at) = at isa Expr && at.head === :tuple ? String.(string.(at.args)) :
                  at === :(()) ? String[] : error("argtypes not a tuple: $at")

function emit(io, sym::Symbol, ret, at)
    ats = argtypes_of(at)
    sig  = join((CTYPE[t] for t in ats), ", ")
    args = join((argexpr(t, k-1) for (k, t) in enumerate(ats)), ", ")
    rs   = string(ret)
    println(io, "  if (sym == \"$sym\") {")
    if rs == "Cvoid"
        println(io, "    ((void(*)($sig))fn)($args);")
        println(io, "    reply_void(out);")
    elseif rs == "Cint"
        println(io, "    int r = ((int(*)($sig))fn)($args);")
        println(io, "    reply_int(out, r);")
    elseif rs == "Cdouble"
        println(io, "    double r = ((double(*)($sig))fn)($args);")
        println(io, "    reply_double(out, r);")
    elseif rs == "Handle"
        println(io, "    void* r = ((void*(*)($sig))fn)($args);")
        println(io, "    reply_ptr(out, r);")
    else
        error("no C return handling for $rs")
    end
    println(io, "    return true;")
    println(io, "  }")
end

function main()
    out = length(ARGS) >= 1 ? ARGS[1] :
        normpath(joinpath(@__DIR__, "..", "src", "enzomodules_worker_dispatch.inc"))
    m = EnzoLib.manifest()
    h = EnzoLib.contract_hash()
    open(out, "w") do io
        println(io, "// AUTO-GENERATED from Vespa.jl/lib/EnzoLib/src/session.jl @xcall sites")
        println(io, "// by EnzoModules/tools/gen_worker_dispatch.jl — DO NOT EDIT BY HAND.")
        println(io, "// Regenerate after changing the bridge surface; the baked contract hash")
        println(io, "// below must match EnzoLib.contract_hash() or the worker handshake fails.")
        println(io, "#define WORKER_CONTRACT_HASH 0x$(string(h; base=16))ULL")
        println(io, "// $(length(m)) bridge symbols")
        println(io)
        println(io, "// dispatch(sym, fn, a, out): fn = dlsym'd C function pointer; a = decoded")
        println(io, "// args (scalars + mmap'd buffer pointers in a[k].p); out = reply buffer.")
        println(io, "static bool worker_dispatch(const std::string& sym, void* fn,")
        println(io, "                            std::vector<Arg>& a, std::string& out) {")
        for k in sort!(collect(keys(m)); by = string)
            ret, at = m[k]
            emit(io, k, ret, at)
        end
        println(io, "  return false;  // unknown symbol")
        println(io, "}")
    end
    println("wrote $out  ($(length(m)) symbols, contract 0x$(string(h; base=16)))")
end

main()
