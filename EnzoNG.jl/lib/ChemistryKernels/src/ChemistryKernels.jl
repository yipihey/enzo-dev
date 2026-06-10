"""
    ChemistryKernels

A KernelAbstractions.jl port of Enzo's non-equilibrium primordial chemistry and
radiative-cooling network -- the Fortran chain

    calc_rates.F  →  cool1d_multi.F  →  solve_rate_cool.F

(Abel, Anninos, Zhang & Norman 1997; Anninos et al. 1997) -- written **once** and
run on both the CPU (the parity oracle) and the Metal GPU.

Design contract (identical to `PPMKernels`):

  * **One source, three networks, two devices.** Every compute kernel is a
    precision-generic `@kernel` parameterised on the element type `T` *and* on
    the network size through `Val{NSP}` (`NSP ∈ (6, 9, 12)`):

        6  -> de HI HII HeI HeII HeIII                  (MultiSpecies = 1)
        9  -> + HM H2I H2II                             (MultiSpecies = 2)
        12 -> + DI DII HDI                              (MultiSpecies = 3)

    The CPU backend runs it in `Float64` (to certify the algorithm *is* Enzo's,
    bit-tight against the Fortran golden fixtures) and in `Float32`; the Metal
    backend runs it only in `Float32`. f32 CPU↔Metal agreement is the parity gate.
    Stiff chemistry genuinely wants f64 -- f32 on Metal is the same accuracy
    compromise GPU-Grackle makes, and the certification layer flags where it bites.

  * **Backend by name.** `backend(:cpu)` always works; `backend(:metal)` resolves
    only after `using Metal` has loaded the package extension (Apple Silicon).
    Allocation/host-transfer go through [`device_zeros`](@ref) / [`to_device`](@ref)
    / [`to_host`](@ref), which the Metal extension specialises -- the kernels and
    tests never name a concrete array type.

  * **Species-of-arrays layout.** Each species' comoving mass density (Enzo's
    `HI`, `HII`, ... fields, in code density units) is its own contiguous array,
    so the per-cell kernels read/write coalesced columns. A whole grid is a
    `NamedTuple` of these arrays; see [`Species`](@ref).

The network integrator is per-cell *embarrassingly parallel* (one work-item owns
one cell's entire sub-cycled stiff solve), which is why this is the single most
natural GPU target in Enzo -- exactly the structure `pgas2d_dual!` uses for its
row-sequential sweep, here one cell per thread.
"""
module ChemistryKernels

using KernelAbstractions
using Adapt
const KA = KernelAbstractions

export backend, has_backend, device_zeros, to_device, to_host
export RateTables, build_rate_tables, ChemistryUnits
export Species, species_names, nspecies_for
export cool1d_multi!, solve_rate_cool!, edot_cell, edot_scalar

# ── backend registry (verbatim PPMKernels convention) ────────────────────────
const _BACKENDS = Dict{Symbol,Any}(:cpu => CPU())

"Register a KernelAbstractions backend under `name` (used by the Metal extension)."
register_backend!(name::Symbol, be) = (_BACKENDS[name] = be)

"True when backend `name` is available (`:metal` needs `using Metal` first)."
has_backend(name::Symbol) = haskey(_BACKENDS, name)

"""
    backend(name::Symbol = :cpu)

The KernelAbstractions backend registered under `name`. `:cpu` is always
available; `:metal` requires `using Metal` (Apple Silicon) to have loaded the
`ChemistryKernelsMetalExt` extension.
"""
function backend(name::Symbol = :cpu)
    return get(_BACKENDS, name) do
        error("Chemistry backend :$name is not available. " *
              (name === :metal ? "Run `using Metal` first (Apple Silicon only)." :
               "Known backends: $(collect(keys(_BACKENDS)))."))
    end
end

# ── device array helpers (specialised by the Metal extension) ────────────────
"""
    device_zeros(be, T, dims::Dims) -> AbstractArray{T}

A zero-filled array of element type `T` and shape `dims` living on backend `be`.
The CPU default is `zeros(T, dims)`; the Metal extension returns an `MtlArray`.
"""
device_zeros(::CPU, ::Type{T}, dims::Dims) where {T} = zeros(T, dims)

"""
    to_device(be, a, T = eltype(a)) -> AbstractArray{T}

Copy host array `a` onto backend `be`, converting to element type `T`.
"""
function to_device(be, a::AbstractArray, ::Type{T} = eltype(a)) where {T}
    d = device_zeros(be, T, size(a))
    copyto!(d, convert(Array{T}, a))
    return d
end

"`to_host(a)` — a plain host `Array` copy of a device array; synchronizes first."
function to_host(a::AbstractArray)
    KA.synchronize(KA.get_backend(a))
    return Array(a)
end

include("units.jl")        # ChemistryUnits: code↔cgs factors (units.F surrogate)
include("interpolate.jl")  # log-T table index + linear interp (interpolate.F)
include("rates.jl")        # calc_rates.F: k-rates + cooling coefficient tables
include("cooling.jl")      # cool1d_multi.F: per-cell edot assembly
include("network.jl")      # solve_rate_cool.F: sub-cycled semi-implicit BDF step
include("generic.jl")      # arbitrary mass-action network executor (KROME target)
include("adapt.jl")        # Adapt rules so the SoA structs reach device kernels

end # module
