# ── Sub-cycled semi-implicit network step (port of solve_rate_cool.F) ─────────
# One work-item owns ONE cell and runs that cell's entire sub-cycled stiff solve
# — the per-cell-independent structure that makes chemistry the most natural GPU
# target in Enzo (cf. `pgas2d_dual!`, where one thread owns a whole row).
#
# Each sub-step uses Anninos et al. (1997)'s semi-implicit backward-difference
# update, species by species:
#
#     n_new = (creation·dt + n_old) / (1 + destruction·dt)        (always ≥ 0)
#
# with the creation/destruction (`scoef`/`acoef`) split taken verbatim from
# solve_rate_cool.F's `step_rate`. Electrons are closed by exact charge
# conservation (Enzo's electron equation is built to conserve charge, so this is
# bit-equivalent and unconditionally non-negative). The energy advances
# explicitly, `ge += edot/ρ·dt`, and the adaptive sub-step `dtit` is the
# rate_timestep control: 10% of the electron and internal-energy change, capped
# at the remaining `dt` and `0.5·dt`.

export PhotoRates

"""
    PhotoRates{A}

Optional per-cell radiation coupling fed in from `RadiationKernels` (or a UV
background): photo-ionization rate arrays `kphHI`, `kphHeI`, `kphHeII`, `kdissH2`
and the photo-heating array `photogamma` (energy code units, added to `edot`).
Pass `nothing` for the pure collisional network. Arrays share the grid layout.
"""
struct PhotoRates{A}
    kphHI::A; kphHeI::A; kphHeII::A; kdissH2::A; photogamma::A
end

@inline _bdf(scoef::T, acoef::T, nold::T, dt::T) where {T} =
    (scoef * dt + nold) / (one(T) + acoef * dt)

@inline function _electron_density(HII, HeII, HeIII, HM, H2II, ::Val{NSP}) where {NSP}
    T = typeof(HII)
    de = HII + T(0.25) * HeII + T(0.5) * HeIII
    if NSP ≥ 9
        de += T(0.5) * H2II - HM
    end
    return max(de, T(RATE_TINY))
end

# The whole-cell solve, as a plain device function so both CPU and Metal inline it.
@inline function _solve_cell!(sp, idx, dens::T, ge_in::T, dt::T, gamma::T,
                              rt, units, temperature_units::T,
                              comp1::T, comp2::T, gammah::T, photo, has_photo::Bool,
                              conserve::Bool, itmax::Int, ::Val{NSP}) where {T,NSP}
    @inbounds begin
        de  = sp.de[idx]
        HI  = sp.HI[idx];  HII  = sp.HII[idx]
        HeI = sp.HeI[idx]; HeII = sp.HeII[idx]; HeIII = sp.HeIII[idx]
        HM   = NSP ≥ 9 ? sp.HM[idx]   : zero(T)
        H2I  = NSP ≥ 9 ? sp.H2I[idx]  : zero(T)
        H2II = NSP ≥ 9 ? sp.H2II[idx] : zero(T)
        DI   = NSP ≥ 12 ? sp.DI[idx]  : zero(T)
        DII  = NSP ≥ 12 ? sp.DII[idx] : zero(T)
        HDI  = NSP ≥ 12 ? sp.HDI[idx] : zero(T)
        ge = ge_in

        # conserved nuclei budgets (chemistry only moves atoms between ionization
        # / molecular states), enforced by make_consistent below — Grackle's step.
        Htot  = HI + HII + (NSP ≥ 9 ? HM + H2I + H2II : zero(T))
        Hetot = HeI + HeII + HeIII

        kHI   = has_photo ? photo.kphHI[idx]   : zero(T)
        kHeI  = has_photo ? photo.kphHeI[idx]  : zero(T)
        kHeII = has_photo ? photo.kphHeII[idx] : zero(T)
        kdiss = has_photo ? photo.kdissH2[idx] : zero(T)
        pheat = has_photo ? photo.photogamma[idx] : zero(T)

        ttot = zero(T); q = T(0.25)
        for iter in 1:itmax
            # temperature + table lookup for this sub-step
            ntot = de + HI + HII + (HeI + HeII + HeIII) * q
            if NSP ≥ 9
                ntot += HM + (H2I + H2II) * T(0.5)
            end
            mu = dens / max(ntot, T(RATE_TINY))
            tgas = max((gamma - one(T)) * ge * mu * temperature_units, T(1))
            ii, tdef = table_index(log(tgas), rt.grid)

            k1 = interp(rt.k1, ii, tdef); k2 = interp(rt.k2, ii, tdef)
            k3 = interp(rt.k3, ii, tdef); k4 = interp(rt.k4, ii, tdef)
            k5 = interp(rt.k5, ii, tdef); k6 = interp(rt.k6, ii, tdef)

            # explicit electron rate of change for the sub-step control
            dedot = k1*HI*de - k2*HII*de +
                    (k3*HeI*de - k4*HeII*de + k5*HeII*de - k6*HeIII*de) * q +
                    (has_photo ? kHI*HI + (kHeI*HeI + kHeII*HeII)*q : zero(T))
            # edot from the CURRENT local densities (NOT the stale grid arrays),
            # so cooling tracks ionization within the sub-cycle.
            edot = edot_scalar(de, HI, HII, HeI, HeII, HeIII,
                               (NSP ≥ 9 ? H2I : zero(T)), dens, tgas, rt, ii, tdef,
                               units, comp1, comp2, gammah, Val(NSP)) + pheat

            # adaptive sub-step: 10% of electron & energy change, capped (rate_timestep)
            dtde = abs(dedot) > T(RATE_TINY) ? T(0.1) * de / abs(dedot) : dt
            dte  = abs(edot)  > T(RATE_TINY) ? T(0.1) * abs(dens * ge / edot) : dt
            dtit = min(dtde, dte, dt - ttot, T(0.5) * dt)
            dtit = max(dtit, T(1e-3) * (dt - ttot) + T(RATE_TINY))

            # ── semi-implicit species update (step_rate ordering) ──────────────
            # H ionization balance
            sHI = k2 * HII * de
            aHI = k1 * de + (has_photo ? kHI : zero(T))
            HIn = _bdf(sHI, aHI, HI, dtit)
            sHII = k1 * HIn * de + (has_photo ? kHI * HIn : zero(T))
            aHII = k2 * de
            HIIn = _bdf(sHII, aHII, HII, dtit)
            # He ionization chain
            sHeI = k4 * HeII * de
            aHeI = k3 * de + (has_photo ? kHeI : zero(T))
            HeIn = _bdf(sHeI, aHeI, HeI, dtit)
            sHeII = k3 * HeIn * de + k6 * HeIII * de +
                    (has_photo ? kHeI * HeIn : zero(T))
            aHeII = (k4 + k5) * de + (has_photo ? kHeII : zero(T))
            HeIIn = _bdf(sHeII, aHeII, HeII, dtit)
            sHeIII = k5 * HeIIn * de + (has_photo ? kHeII * HeIIn : zero(T))
            aHeIII = k6 * de
            HeIIIn = _bdf(sHeIII, aHeIII, HeIII, dtit)

            HMn = HM; H2In = H2I; H2IIn = H2II
            if NSP ≥ 9
                k7  = interp(rt.k7, ii, tdef);  k8  = interp(rt.k8, ii, tdef)
                k9  = interp(rt.k9, ii, tdef);  k10 = interp(rt.k10, ii, tdef)
                k11 = interp(rt.k11, ii, tdef); k12 = interp(rt.k12, ii, tdef)
                k13 = interp(rt.k13, ii, tdef); k14 = interp(rt.k14, ii, tdef)
                k15 = interp(rt.k15, ii, tdef); k16 = interp(rt.k16, ii, tdef)
                k17 = interp(rt.k17, ii, tdef); k18 = interp(rt.k18, ii, tdef)
                k19 = interp(rt.k19, ii, tdef)
                # H-  : form k7·HI·de ; destroy k8·HI + k15·HI + k16·HII + k17·HII + k14·de
                sHM = k7 * HIn * de
                aHM = (k8 + k15) * HIn + (k16 + k17) * HIIn + k14 * de
                HMn = _bdf(sHM, aHM, HM, dtit)
                # H2+ : form k9·HI·HII + k11·H2I·HII + k17·HM·HII ; destroy k10·HI + k18·de + k19·HM
                sH2II = k9 * HIn * HIIn + k11 * H2I * HIIn + k17 * HMn * HIIn
                aH2II = k10 * HIn + k18 * de + k19 * HMn
                H2IIn = _bdf(sH2II, aH2II, H2II, dtit)
                # H2  : form k8·HM·HI + k10·H2II·HI + k19·H2II·HM ; destroy k11·HII + k12·de + k13·HI
                # (+ Lyman-Werner photodissociation kdiss from the radiation coupling)
                sH2I = k8 * HMn * HIn + k10 * H2IIn * HIn + k19 * H2IIn * HMn
                aH2I = k11 * HIIn + k12 * de + k13 * HIn + (has_photo ? kdiss : zero(T))
                H2In = _bdf(sH2I, aH2I, H2I, dtit)
            end

            # ── deuterium / HD network (NSP=12; Enzo solve_rate_cool.F block D) ─
            DIn = DI; DIIn = DII; HDIn = HDI
            if NSP ≥ 12
                k50 = interp(rt.k50, ii, tdef); k51 = interp(rt.k51, ii, tdef)
                k52 = interp(rt.k52, ii, tdef); k53 = interp(rt.k53, ii, tdef)
                k54 = interp(rt.k54, ii, tdef); k55 = interp(rt.k55, ii, tdef)
                k56 = interp(rt.k56, ii, tdef)
                # D shares H's collisional ion/recomb (k1/k2) and photo-ionization.
                # DI : form k2·DII·de + k51·DII·HI + 2·k55·HDI·HI/3 ;
                #      destroy k1·de + k50·HII + k54·H2I/2 + k56·HM (+ photo)
                sDI = k2*DIIn*de + k51*DIIn*HIn + T(2)*k55*HDI*HIn/T(3)
                aDI = k1*de + k50*HIIn + k54*H2In/T(2) + k56*HMn +
                      (has_photo ? kHI : zero(T))
                DIn = _bdf(sDI, aDI, DI, dtit)
                # DII: form k1·DI·de + k50·HII·DI + 2·k53·HII·HDI/3 (+ photo·DI) ;
                #      destroy k2·de + k51·HI + k52·H2I/2
                sDII = k1*DIn*de + k50*HIIn*DIn + T(2)*k53*HIIn*HDI/T(3) +
                       (has_photo ? kHI*DIn : zero(T))
                aDII = k2*de + k51*HIn + k52*H2In/T(2)
                DIIn = _bdf(sDII, aDII, DII, dtit)
                # HDI: form 3·(k52·DII·H2I/4 + k54·DI·H2I/4 + 2·k56·DI·HM/2) ;
                #      destroy k53·HII + k55·HI
                sHDI = T(3)*(k52*DIIn*H2In/T(4) + k54*DIn*H2In/T(4) +
                             T(2)*k56*DIn*HMn/T(2))
                aHDI = k53*HIIn + k55*HIn
                HDIn = _bdf(sHDI, aHDI, HDI, dtit)
            end

            # commit, close electrons by charge conservation, advance energy/clock
            HI = max(HIn, T(RATE_TINY));  HII = max(HIIn, T(RATE_TINY))
            HeI = max(HeIn, T(RATE_TINY)); HeII = max(HeIIn, T(RATE_TINY))
            HeIII = max(HeIIIn, T(RATE_TINY))
            HM = max(HMn, T(RATE_TINY)); H2I = max(H2In, T(RATE_TINY))
            H2II = max(H2IIn, T(RATE_TINY))
            DI = max(DIn, T(RATE_TINY)); DII = max(DIIn, T(RATE_TINY))
            HDI = max(HDIn, T(RATE_TINY))
            de = _electron_density(HII, HeII, HeIII, HM, H2II, Val(NSP))
            ge = max(ge + edot / dens * dtit, T(RATE_TINY))

            ttot += dtit
            ttot ≥ dt * (one(T) - T(1e-6)) && break
        end

        # make_consistent (optional): renormalize H/He species to their fixed
        # nuclei budgets so the independent semi-implicit updates conserve exactly.
        # Enzo's native solve_rate_cool does NOT do this -- pass `conserve=false`
        # for bit-exact Fortran parity, `true` (default) for robustness.
        if conserve
            Hsum = HI + HII + (NSP ≥ 9 ? HM + H2I + H2II : zero(T))
            sH = Htot / max(Hsum, T(RATE_TINY))
            HI *= sH; HII *= sH
            if NSP ≥ 9
                HM *= sH; H2I *= sH; H2II *= sH
            end
            Hesum = HeI + HeII + HeIII
            sHe = Hetot / max(Hesum, T(RATE_TINY))
            HeI *= sHe; HeII *= sHe; HeIII *= sHe
            de = _electron_density(HII, HeII, HeIII, HM, H2II, Val(NSP))
        end

        # write the updated state back
        sp.de[idx] = de
        sp.HI[idx] = HI;   sp.HII[idx] = HII
        sp.HeI[idx] = HeI; sp.HeII[idx] = HeII; sp.HeIII[idx] = HeIII
        if NSP ≥ 9
            sp.HM[idx] = HM; sp.H2I[idx] = H2I; sp.H2II[idx] = H2II
        end
        if NSP ≥ 12
            sp.DI[idx] = DI; sp.DII[idx] = DII; sp.HDI[idx] = HDI
        end
        return ge
    end
end

@kernel function _solve_rate_cool_kernel!(ge, @Const(dens), sp, rt, units,
                                          gamma, temperature_units, comp1, comp2,
                                          gammah, photo, has_photo::Bool, conserve::Bool,
                                          dt, itmax::Int, ::Val{NSP},
                                          i1::Int, j1::Int, idim::Int) where {NSP}
    gi, gj = @index(Global, NTuple)
    i = i1 + gi - 1; j = j1 + gj - 1
    idx = (j - 1) * idim + i
    T = eltype(ge)
    @inbounds ge[idx] = _solve_cell!(sp, idx, dens[idx], ge[idx], T(dt), T(gamma),
                                     rt, units, T(temperature_units),
                                     T(comp1), T(comp2), T(gammah), photo, has_photo,
                                     conserve, itmax, Val(NSP))
end

"""
    solve_rate_cool!(ge, dens, sp::Species, rt::RateTables, units;
                     nspecies, gamma, temperature_units, dt,
                     idim, i1, i2, j1=1, j2=1, itmax=10000,
                     comp1=0, comp2=0, gammah=0, photo=nothing) -> ge

Advance the species densities in `sp` and the specific internal energy `ge` by one
hydro step `dt` (code units), sub-cycling each cell's stiff chemistry+cooling
independently (`solve_rate_cool.F`). `nspecies ∈ (6,9,12)` selects the network;
pass a [`PhotoRates`](@ref) for `photo` to couple in radiative-transfer
ionization/heating, or `nothing` for the collisional-only network. Updates `sp`
in place and returns the updated `ge`.
"""
function solve_rate_cool!(ge, dens, sp::Species, rt, units;
                          nspecies::Integer, gamma::Real, temperature_units::Real,
                          dt::Real, idim::Integer, i1::Integer, i2::Integer,
                          j1::Integer = 1, j2::Integer = 1, itmax::Integer = 10000,
                          comp1::Real = 0, comp2::Real = 0, gammah::Real = 0,
                          photo = nothing, conserve::Bool = true)
    be = KA.get_backend(ge)
    ni = Int(i2 - i1 + 1); nj = Int(j2 - j1 + 1)
    has_photo = photo !== nothing
    ph = has_photo ? photo : PhotoRates(ge, ge, ge, ge, ge)  # dummy, never read
    _solve_rate_cool_kernel!(be)(ge, dens, sp, rt, units, gamma, temperature_units,
                                 comp1, comp2, gammah, ph, has_photo, conserve, dt,
                                 Int(itmax), Val(Int(nspecies)),
                                 Int(i1), Int(j1), Int(idim); ndrange = (ni, nj))
    return ge
end
