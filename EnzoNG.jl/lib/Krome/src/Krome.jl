"""
    Krome

A Julia front-end for the **KROME** astrochemistry package (Grassi et al. 2014;
<https://bitbucket.org/tgrassi/krome>, GPL-3.0). It reads KROME `react_*` network
files, compiles their Fortran rate expressions into Julia, and lowers the network
onto [`ChemistryKernels.GenericNetwork`](@ref) — so any KROME network runs on the
CPU and the GPU through KernelAbstractions, the same single-source contract as
EnzoNG's hand-written primordial network.

Pipeline:

    react_* file ──parse_krome──▶ KromeNetwork ──to_generic──▶ GenericNetwork
                                   (species + reactions,        (tabulated rates +
                                    Tgas→k closures)             stoichiometry, CPU/GPU)

The upstream KROME tree is vendored as a git submodule at `extern/krome` (under
its own GPL-3.0 license); [`krome_network_path`](@ref) resolves its bundled
`networks/react_*` databases.

Units: KROME rate coefficients are cgs (cm³/s for two-body). The generic executor
therefore expects **number densities (cm⁻³)** and **seconds**, not Enzo code units
— distinct from `ChemistryKernels.solve_rate_cool!`, which works in code units.
"""
module Krome

using ChemistryKernels
using KernelAbstractions

export KromeNetwork, KromeReaction, parse_krome, to_generic, direct_step!
export krome_network_path, list_krome_networks, compile_rate

include("expr.jl")    # Fortran rate-expression → Julia closure
include("parse.jl")   # react_* file → KromeNetwork IR
include("build.jl")   # KromeNetwork → ChemistryKernels.GenericNetwork

"""
    krome_network_path(name) -> String

Absolute path to the vendored KROME network file `name` (e.g. `"react_primordial"`)
under `extern/krome/networks`. Errors with a hint if the submodule is not checked
out (`git submodule update --init`).
"""
function krome_network_path(name::AbstractString)
    root = normpath(joinpath(@__DIR__, "..", "extern", "krome", "networks"))
    isdir(root) || error("KROME submodule not initialised — run " *
        "`git submodule update --init EnzoNG.jl/lib/Krome/extern/krome` (path: $root).")
    p = joinpath(root, name)
    isfile(p) || error("KROME network '$name' not found under $root.")
    return p
end

"`list_krome_networks()` — names of the vendored `react_*` files (submodule)."
function list_krome_networks()
    root = normpath(joinpath(@__DIR__, "..", "extern", "krome", "networks"))
    isdir(root) || return String[]
    return sort(filter(f -> startswith(f, "react_"), readdir(root)))
end

end # module
