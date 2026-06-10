"""
    Grackle

A Julia front-end for the **Grackle** chemistry & cooling library
(<https://github.com/grackle-project/grackle>, Enzo Public License). Grackle
descends from Enzo's `solve_rate_cool.F`/`cool1d_multi.F` — the same primordial
network already ported in `ChemistryKernels` — so the new content wired here is
Grackle's **Cloudy metal-line cooling**:

  * [`load_cloudy_table`](@ref) reads Grackle's Cloudy HDF5 tables into a
    `ChemistryKernels.CloudyTable`, which the KA `cloudy_metal_cooling!` kernel
    interpolates on the CPU and the GPU (the hot path).
  * [`libgrackle`](@ref) / [`grackle_version`](@ref) bind the built `libgrackle.so`
    as the ground-truth oracle for parity/benchmarks.

Upstream Grackle is vendored as a git submodule at `extern/grackle`; build it with
`extern/grackle` + CMake (see `deps/build_grackle.sh`) and point `Grackle.jl` at
the resulting library via the `GRACKLE_LIB` environment variable or
[`set_libgrackle!`](@ref).
"""
module Grackle

using HDF5
using ChemistryKernels
using KernelAbstractions

export load_cloudy_table, libgrackle, grackle_version, set_libgrackle!, grackle_available

# ── libgrackle discovery ──────────────────────────────────────────────────────
const _LIB = Ref{String}(get(ENV, "GRACKLE_LIB", ""))
_libpath::String = ""     # ccall library must be a (typed) global, not a local

"Point Grackle.jl at a built `libgrackle.so`."
set_libgrackle!(path::AbstractString) = (_LIB[] = String(path); nothing)

function _find_lib()
    isempty(_LIB[]) || return _LIB[]
    for cand in (joinpath(@__DIR__, "..", "deps", "build", "libgrackle.so"),
                 get(ENV, "GRACKLE_LIB", ""))
        (!isempty(cand) && isfile(cand)) && return (_LIB[] = cand)
    end
    return ""
end

"`libgrackle()` — path to the built shared library, or `\"\"` if not found."
libgrackle() = _find_lib()

"True when a built `libgrackle.so` is locatable (binding/oracle available)."
grackle_available() = !isempty(_find_lib())

"""
    grackle_version() -> NamedTuple

Query the built `libgrackle`'s `get_grackle_version()` (version / branch / revision).
Errors if the library is not found ([`grackle_available`](@ref) is false).
"""
function grackle_version()
    isempty(_find_lib()) && error("libgrackle not found — build it and set " *
                          "GRACKLE_LIB (see deps/build_grackle.sh).")
    global _libpath = _LIB[]          # ccall lib must be a global, not a local
    # get_grackle_version returns a struct of 3 char* (version, branch, revision)
    gv = ccall((:get_grackle_version, _libpath), NTuple{3,Cstring}, ())
    return (; version = unsafe_string(gv[1]), branch = unsafe_string(gv[2]),
            revision = unsafe_string(gv[3]))
end

include("cloudy_io.jl")   # HDF5 Cloudy-table loader → ChemistryKernels.CloudyTable

end # module
