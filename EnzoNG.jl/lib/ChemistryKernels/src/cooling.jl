# ── Per-cell cooling/heating rate (port of cool1d_multi.F) ────────────────────
# `edot` is the net volumetric energy-change rate (heating − cooling) in cooling
# code units. The terms are assembled exactly as cool1d_multi.F's big sum:
# collisional excitation + ionization + recombination + bremsstrahlung + Compton
# (+ H2/HD line cooling for the 9/12-species networks + photoelectric dust
# heating + photo-heating). All coefficients are interpolated from the log(T)
# tables; the He species carry Enzo's `/4` (mass→ number) convention verbatim.
#
# The temperature is recovered from the internal energy and the instantaneous
# mean molecular weight (the sum over species), matching cool1d_multi's `tgas`.

"""
    number_density(sp, idx, ::Val{NSP}) -> n

Total particle number density (in code-density units, i.e. before `·dom`) summed
over the active species at linear cell `idx`: `Σ n_s` with He divided by 4 and H2
counted once. Used for the mean molecular weight / temperature.
"""
@inline function number_density(sp, idx, ::Val{NSP}) where {NSP}
    T = eltype(sp.HI)
    @inbounds begin
        n = sp.de[idx] + sp.HI[idx] + sp.HII[idx] +
            (sp.HeI[idx] + sp.HeII[idx] + sp.HeIII[idx]) * T(0.25)
        if NSP ≥ 9
            n += sp.HM[idx] + (sp.H2I[idx] + sp.H2II[idx]) * T(0.5)
        end
        if NSP ≥ 12
            n += (sp.DI[idx] + sp.DII[idx]) * T(0.5) + sp.HDI[idx] * T(1/3)
        end
    end
    return n
end

"""
    gas_temperature(ge, dens, sp, idx, gamma, units, ::Val{NSP}) -> tgas

Gas temperature (K) from the specific internal energy `ge` (code units) and the
species mix, `T = (γ-1)·μ·m_H/k_B · ge` collapsed into cool1d_multi's form
`tgas = (γ-1)·ge·dens / n_tot · dom_temperature`. Here `ge` is energy per unit
mass; `n_tot` the number density above.
"""
@inline function gas_temperature(ge::T, dens::T, sp, idx, gamma::T,
                                 temperature_units::T, ::Val{NSP}) where {T,NSP}
    ntot = number_density(sp, idx, Val(NSP))
    # μ = dens / ntot ; tgas = (γ-1)·ge·μ·temperature_units
    mu = dens / max(ntot, T(RATE_TINY))
    tgas = (gamma - one(T)) * ge * mu * temperature_units
    return max(tgas, T(1))   # Enzo floors at the table start
end

"""
    edot_cell(sp, idx, dens, tgas, rt, i, tdef, units, comp1, comp2, gammah, ::Val{NSP}) -> edot

The net heating−cooling rate for one cell (cooling code units), assembled term by
term from cool1d_multi.F. `i,tdef` are the table index/weight for `tgas`; `comp1,
comp2` the redshift-scaled Compton coefficient/temperature; `gammah` toggles
photoelectric dust heating. Returns `edot`; the integrator advances
`ge += edot/dens·dt`.
"""
@inline function edot_cell(sp, idx, dens::T, tgas::T, rt, i::Int, tdef::T,
                           units, comp1::T, comp2::T, gammah::T,
                           ::Val{NSP}) where {T,NSP}
    @inbounds begin
        H2I_ = NSP ≥ 9 ? sp.H2I[idx] : zero(T)
        return edot_scalar(sp.de[idx], sp.HI[idx], sp.HII[idx], sp.HeI[idx],
                           sp.HeII[idx], sp.HeIII[idx], H2I_, dens, tgas, rt, i, tdef,
                           units, comp1, comp2, gammah, Val(NSP))
    end
end

"""
    edot_scalar(de, HI, HII, HeI, HeII, HeIII, H2I, dens, tgas, rt, i, tdef,
                units, comp1, comp2, gammah, ::Val{NSP}) -> edot

Net heating−cooling rate from the *current* species values passed as scalars (not
read from a grid array). The sub-cycled integrator MUST use this with its evolving
local densities — computing `edot` from stale array values freezes `n_e` and lets
the gas over-ionize without cooling (the bug this signature prevents).
"""
@inline function edot_scalar(de_::T, HI_::T, HII_::T, HeI_::T, HeII_::T, HeIII_::T,
                             H2I_::T, dens::T, tgas::T, rt, i::Int, tdef::T,
                             units, comp1::T, comp2::T, gammah::T,
                             ::Val{NSP}) where {T,NSP}
    dom     = units.dom
    dom_inv = units.dom_inv
    @inbounds begin
        ceHI   = interp(rt.ceHI, i, tdef);   ceHeI  = interp(rt.ceHeI, i, tdef)
        ceHeII = interp(rt.ceHeII, i, tdef)
        ciHI   = interp(rt.ciHI, i, tdef);   ciHeI  = interp(rt.ciHeI, i, tdef)
        ciHeII = interp(rt.ciHeII, i, tdef); ciHeIS = interp(rt.ciHeIS, i, tdef)
        reHII  = interp(rt.reHII, i, tdef);  reHeII1 = interp(rt.reHeII1, i, tdef)
        reHeII2 = interp(rt.reHeII2, i, tdef); reHeIII = interp(rt.reHeIII, i, tdef)
        brem   = interp(rt.brem, i, tdef)

        q = T(0.25)   # He /4 convention
        edot =
            # collisional excitation
            - ceHI * HI_ * de_
            - ceHeI * HeII_ * de_^2 * dom * q
            - ceHeII * HeII_ * de_ * q
            # collisional ionization
            - ciHI * HI_ * de_
            - ciHeI * HeI_ * de_ * q
            - ciHeII * HeII_ * de_ * q
            - ciHeIS * HeII_ * de_^2 * dom * q
            # recombination
            - reHII * HII_ * de_
            - reHeII1 * HeII_ * de_ * q
            - reHeII2 * HeII_ * de_ * q
            - reHeIII * HeIII_ * de_ * q
            # Compton (cooling/heating off the CMB) + X-ray Compton
            - comp1 * (tgas - comp2) * de_ * dom_inv
            - rt.comp_xray * (tgas - rt.comp_temp) * de_ * dom_inv
            # bremsstrahlung
            - brem * (HII_ + HeII_ * q + HeIII_) * de_
            # photoelectric dust heating
            + gammah * rt.gammaha * (HI_ + HII_) * dom_inv

        # ── H2 / HD line cooling (9- and 12-species networks) ─────────────────
        if NSP ≥ 9
            galdl = interp(rt.gpldl, i, tdef)              # low-density limit
            gahdl = interp(rt.gphdl, i, tdef)              # high-density (roth)
            # Lepp & Shull bridging: cool = gahdl / (1 + gahdl/(n·galdl))
            nH = HI_ * dom
            gphdl1 = gahdl / max(galdl * nH, T(RATE_TINY))
            h2cool = gahdl / (one(T) + gphdl1)
            edot -= h2cool * H2I_ * q
        end
    end
    return edot
end

# ── kernel: fill an edot field over the active region (diagnostic / Grackle-style)
@kernel function _cool1d_multi_kernel!(edot, @Const(ge), @Const(dens),
                                       sp, rt, units, gamma, temperature_units,
                                       comp1, comp2, gammah, ::Val{NSP},
                                       i1::Int, j1::Int, idim::Int) where {NSP}
    gi, gj = @index(Global, NTuple)
    i = i1 + gi - 1; j = j1 + gj - 1
    idx = (j - 1) * idim + i
    T = eltype(edot)
    @inbounds begin
        d = dens[idx]
        tgas = gas_temperature(ge[idx], d, sp, idx, T(gamma), T(temperature_units), Val(NSP))
        ii, tdef = table_index(log(tgas), rt.grid)
        edot[idx] = edot_cell(sp, idx, d, tgas, rt, ii, tdef, units,
                              T(comp1), T(comp2), T(gammah), Val(NSP))
    end
end

"""
    cool1d_multi!(edot, ge, dens, sp::Species, rt::RateTables, units;
                  nspecies, gamma, temperature_units, idim, i1, i2, j1=1, j2=1,
                  comp1=0, comp2=0, gammah=0) -> edot

Fill `edot` with the per-cell net heating−cooling rate over the active region
(`cool1d_multi.F`). All arrays live on one backend; `nspecies ∈ (6,9,12)` selects
the network. This is the standalone cooling-rate evaluator (used for cooling-time
diagnostics and parity tests); [`solve_rate_cool!`](@ref) calls the same
`edot_cell` inside its sub-cycle.
"""
function cool1d_multi!(edot, ge, dens, sp::Species, rt, units;
                       nspecies::Integer, gamma::Real, temperature_units::Real,
                       idim::Integer, i1::Integer, i2::Integer,
                       j1::Integer = 1, j2::Integer = 1,
                       comp1::Real = 0, comp2::Real = 0, gammah::Real = 0)
    be = KA.get_backend(edot)
    ni = Int(i2 - i1 + 1); nj = Int(j2 - j1 + 1)
    _cool1d_multi_kernel!(be)(edot, ge, dens, sp, rt, units, gamma, temperature_units,
                              comp1, comp2, gammah, Val(Int(nspecies)),
                              Int(i1), Int(j1), Int(idim); ndrange = (ni, nj))
    return edot
end
