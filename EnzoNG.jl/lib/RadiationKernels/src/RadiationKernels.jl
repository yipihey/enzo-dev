"""
    RadiationKernels

A KernelAbstractions.jl two-moment (M1) radiative-transfer solver for EnzoNG --
the GPU-native counterpart to Enzo's adaptive ray tracing -- written **once** and
run on the CPU (the parity oracle) and the Metal GPU.

The state per photon group `g` is the photon number density `N[g]` and the photon
number-flux vector `F[g] = (Fx, Fy, Fz)`. The system

    ∂ₜ N   + ∇·F        =  Ṅ_source − c Σ_s n_s σ_{s,g} N
    ∂ₜ F   + c² ∇·P     =        − c Σ_s n_s σ_{s,g} F

is closed by the Levermore **M1** Eddington tensor [`eddington_factor`](@ref)
(Aubert & Teyssier 2008; Rosdahl et al. 2013 / RAMSES-RT), transported by a
dimensionally-split global-Lax-Friedrichs (GLF) scheme [`m1_sweep!`](@ref), and
coupled to the chemistry network through per-cell photo-ionization/heating
[`photo_rates!`](@ref). A reduced speed of light `c_red` is supported.

Design contract matches `PPMKernels`/`ChemistryKernels`: precision-generic
`@kernel`s, a backend registry (`:cpu` always, `:metal` via the package
extension), and device array helpers the Metal extension specialises. Group count
is a call-site `Val{NG}`, so 1-, 3-, … group transport is one kernel.
"""
module RadiationKernels

using KernelAbstractions
using Adapt
const KA = KernelAbstractions

export backend, has_backend, device_zeros, to_device, to_host
export eddington_factor, RadField, m1_sweep!, m1_step!
export PhotoCrossSections, photo_rates!, SPEED_OF_LIGHT

const SPEED_OF_LIGHT = 2.99792458e10   # cm/s

# ── backend registry (verbatim PPMKernels convention) ────────────────────────
const _BACKENDS = Dict{Symbol,Any}(:cpu => CPU())
register_backend!(name::Symbol, be) = (_BACKENDS[name] = be)
has_backend(name::Symbol) = haskey(_BACKENDS, name)

"""
    backend(name::Symbol = :cpu)

The KernelAbstractions backend registered under `name`. `:cpu` always available;
`:metal` requires `using Metal` (loads `RadiationKernelsMetalExt`).
"""
function backend(name::Symbol = :cpu)
    return get(_BACKENDS, name) do
        error("Radiation backend :$name is not available. " *
              (name === :metal ? "Run `using Metal` first (Apple Silicon only)." :
               "Known backends: $(collect(keys(_BACKENDS)))."))
    end
end

device_zeros(::CPU, ::Type{T}, dims::Dims) where {T} = zeros(T, dims)
function to_device(be, a::AbstractArray, ::Type{T} = eltype(a)) where {T}
    d = device_zeros(be, T, size(a)); copyto!(d, convert(Array{T}, a)); return d
end
function to_host(a::AbstractArray)
    KA.synchronize(KA.get_backend(a)); return Array(a)
end

include("closure.jl")     # M1 Eddington-tensor closure
include("transport.jl")   # GLF dimensionally-split moment transport
include("coupling.jl")    # per-cell photo-ionization/heating → chemistry

# ── Adapt rules so RadField reaches device kernels (identity on CPU) ──────────
Adapt.adapt_structure(to, f::RadField) = RadField(
    Adapt.adapt(to, f.N), Adapt.adapt(to, f.Fx),
    Adapt.adapt(to, f.Fy), Adapt.adapt(to, f.Fz))

end # module
