# ── Units + species layout ────────────────────────────────────────────────────
# Mirrors the unit bookkeeping Enzo threads through calc_rates / cool1d_multi /
# solve_rate_cool. Enzo stores the chemistry fields as *comoving mass densities*
# in code units; rate equations need *proper number densities*, hence the `dom`
# (density→ number/cm³ of H) and `coolunit` factors carried here.

"Physical constants in cgs (match Enzo's phys_const.def to the digits used there)."
module PhysConst
const mass_h  = 1.67262171e-24   # g  (proton mass, Enzo's mh)
const mass_e  = 9.10938215e-28   # g
const kboltz  = 1.3806504e-16    # erg/K
const everg   = 1.60217653e-12   # erg/eV
const tevk    = 1.1604505e4      # K per eV  (everg/kboltz)
const pi_val  = 3.14159265358979
end

"""
    ChemistryUnits{T}

Code↔cgs conversion factors for the chemistry network at a given snapshot, the
exact set `calc_rates`/`cool1d_multi`/`solve_rate_cool` consume. Build with
[`chemistry_units`](@ref) from the simulation's base units.

  * `dom`       — multiply a code *density* by `dom` to get proper H number
                  density per cm³: `dom = urho·aye³ / m_H`.
  * `coolunit`  — cooling-rate code unit: `(aye⁵·xbase²·m_H²)/(tbase³·dens)`.
  * `kunit`     — two-body rate-coefficient code unit `urho·aye³/(m_H·tbase)`.
  * `coolunit_inv`, `dom_inv` — cached reciprocals used in the per-cell hot loop.
  * `time_to_cgs` — `tbase` (code time → s) for the sub-cycle clock.
"""
struct ChemistryUnits{T}
    dom::T
    dom_inv::T
    coolunit::T
    coolunit_inv::T
    kunit::T
    time_to_cgs::T
    aye::T
end

"""
    chemistry_units(T; urho, uxyz, utim, aye) -> ChemistryUnits{T}

Assemble the chemistry unit factors from Enzo's base units (`urho` density,
`uxyz` length, `utim` time, in cgs) at expansion factor `aye`. Follows the
definitions at the top of `calc_rates.F` / `cool1d_multi.F`.
"""
function chemistry_units(::Type{T}; urho, uxyz, utim, aye) where {T}
    mh   = PhysConst.mass_h
    dom  = urho * aye^3 / mh
    # coolunit = (aye^5 * xbase1^2 * mh^2) / (tbase1^3 * dbase1), with
    # dbase1 = urho*aye^3 (proper) and xbase1 = uxyz/aye (comoving→ proper length).
    xbase1 = uxyz / aye
    dbase1 = urho * aye^3
    tbase1 = utim
    coolunit = (aye^5 * xbase1^2 * mh^2) / (tbase1^3 * dbase1)
    kunit    = (urho * aye^3) / (mh * tbase1)
    return ChemistryUnits{T}(T(dom), T(1/dom), T(coolunit), T(1/coolunit),
                             T(kunit), T(tbase1), T(aye))
end

# ── Species layout ────────────────────────────────────────────────────────────
# The field ordering matches Enzo's IdentifySpeciesFields. `de` is the electron
# density carried in *proton-mass* units (Enzo convention: divide He species by 4
# in the rate equations, which the kernels reproduce verbatim).

const SPECIES6  = (:de, :HI, :HII, :HeI, :HeII, :HeIII)
const SPECIES9  = (SPECIES6..., :HM, :H2I, :H2II)
const SPECIES12 = (SPECIES9..., :DI, :DII, :HDI)

"""
    species_names(nsp) -> Tuple{Vararg{Symbol}}

The ordered field names for a `nsp`-species network (`nsp ∈ (6, 9, 12)`).
"""
function species_names(nsp::Integer)
    nsp == 6  && return SPECIES6
    nsp == 9  && return SPECIES9
    nsp == 12 && return SPECIES12
    error("Unsupported network size $nsp (expected 6, 9 or 12).")
end

"`nspecies_for(ms)` — species count for Enzo's `MultiSpecies` flag (1→6, 2→9, 3→12)."
nspecies_for(ms::Integer) = (6, 9, 12)[ms]

"""
    Species(arrays::NamedTuple)

A thin wrapper over the per-species density arrays of one grid. Construct from a
`NamedTuple` whose keys are a prefix of [`SPECIES12`](@ref); `nspecies(s)` returns
how many are present. All arrays must live on the same backend and share shape.
"""
struct Species{NT<:NamedTuple}
    fields::NT
end
Species(; kwargs...) = Species(values(kwargs))

Base.getproperty(s::Species, k::Symbol) =
    k === :fields ? getfield(s, :fields) : getproperty(getfield(s, :fields), k)
nspecies(s::Species) = length(getfield(s, :fields))
Base.keys(s::Species) = keys(getfield(s, :fields))
