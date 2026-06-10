# ── Radiation → chemistry coupling (port of RadiativeTransferIonization.C) ────
# Turns the local radiation moment field (photon number densities per group) into
# the per-cell photo-ionization rates `kphHI/kphHeI/kphHeII`, the H2 dissociation
# rate `kdissH2`, and the photo-heating rate `photogamma` that the
# ChemistryKernels network consumes as a `PhotoRates`. For group `g` with
# group-averaged cross section σ_{s,g} (cm²) and mean excess energy ΔE_{s,g} (the
# photon energy above the ionization threshold of species s):
#
#     kph_s   = c Σ_g σ_{s,g} N_g                       [1/s, per atom]
#     Γ_heat  = c Σ_g Σ_s n_s σ_{s,g} ΔE_{s,g} N_g      [erg/s/cm³ → energy units]
#
# matching RadiativeTransferIonization's per-absorption rate/heating bookkeeping,
# but evaluated locally from the M1 field instead of along rays.

"""
    PhotoCrossSections{T}

Group-averaged photo cross sections and excess energies for one photon group,
in cgs. `sigmaHI/HeI/HeII` are σ (cm²); `eHI/HeI/HeII` the mean photo-electron
energy above threshold (erg); `sigmaH2` the Lyman-Werner H2-dissociation cross
section. Build a `Vector{PhotoCrossSections}` (one per group) for multi-group.
"""
struct PhotoCrossSections{T}
    sigmaHI::T;  eHI::T
    sigmaHeI::T; eHeI::T
    sigmaHeII::T; eHeII::T
    sigmaH2::T
end

@kernel function _photo_rates_kernel!(kphHI, kphHeI, kphHeII, kdissH2, photogamma,
                                      @Const(N), @Const(nHI), @Const(nHeI), @Const(nHeII),
                                      c, σHI, eHI, σHeI, eHeI, σHeII, eHeII, σH2,
                                      energy_units, accumulate::Bool)
    idx = @index(Global, Linear)
    T = eltype(kphHI)
    @inbounds begin
        Ng = N[idx]
        cN = T(c) * Ng
        dkHI   = T(σHI)   * cN
        dkHeI  = T(σHeI)  * cN
        dkHeII = T(σHeII) * cN
        dkH2   = T(σH2)   * cN
        # photo-heating: c · N · Σ_s n_s σ_s ΔE_s, converted to energy code units
        heat = cN * (nHI[idx]  * T(σHI)   * T(eHI) +
                     nHeI[idx] * T(σHeI)  * T(eHeI) +
                     nHeII[idx]* T(σHeII) * T(eHeII)) * T(energy_units)
        if accumulate
            kphHI[idx]     += dkHI
            kphHeI[idx]    += dkHeI
            kphHeII[idx]   += dkHeII
            kdissH2[idx]   += dkH2
            photogamma[idx]+= heat
        else
            kphHI[idx]     = dkHI
            kphHeI[idx]    = dkHeI
            kphHeII[idx]   = dkHeII
            kdissH2[idx]   = dkH2
            photogamma[idx]= heat
        end
    end
end

"""
    photo_rates!(kphHI, kphHeI, kphHeII, kdissH2, photogamma,
                 fields::Vector{RadField}, xs::Vector{PhotoCrossSections},
                 nHI, nHeI, nHeII; c=SPEED_OF_LIGHT, energy_units=1.0,
                 ndrange) -> (kphHI, …)

Fill the per-cell photo-ionization / dissociation rates and photo-heating from the
multi-group M1 radiation `fields` and their `xs` cross sections, given the proper
number densities `nHI/nHeI/nHeII` (for the heating sum). The first group writes,
subsequent groups accumulate, so the five output arrays hold the total over
groups -- ready to wrap in `ChemistryKernels.PhotoRates`. `c` may be the reduced
speed of light to match the transport.
"""
function photo_rates!(kphHI, kphHeI, kphHeII, kdissH2, photogamma,
                      fields::Vector{<:RadField}, xs::Vector{<:PhotoCrossSections},
                      nHI, nHeI, nHeII;
                      c::Real = SPEED_OF_LIGHT, energy_units::Real = 1.0)
    @assert length(fields) == length(xs) "one cross-section set per photon group"
    be = KA.get_backend(kphHI)
    nd = length(kphHI)
    for (g, (fld, x)) in enumerate(zip(fields, xs))
        _photo_rates_kernel!(be)(kphHI, kphHeI, kphHeII, kdissH2, photogamma,
                                 fld.N, nHI, nHeI, nHeII, c,
                                 x.sigmaHI, x.eHI, x.sigmaHeI, x.eHeI,
                                 x.sigmaHeII, x.eHeII, x.sigmaH2, energy_units,
                                 g > 1; ndrange = nd)
    end
    return kphHI, kphHeI, kphHeII, kdissH2, photogamma
end
