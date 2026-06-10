# ChemistryKernels — CPU correctness/conservation suite.
#
# These are the f64-on-CPU sanity + conservation gates that run headless on any
# CI (no GPU). They are NOT yet the bit-tight Fortran-fixture parity layer: that
# needs `enzomodules_solve_rate_cool` golden outputs (the EnzoLib RPC harness)
# wired in the same way PPMKernels/test certifies against the Fortran. The
# structure here matches that harness so the fixtures drop in unchanged.

using Test
using ChemistryKernels
const CK = ChemistryKernels

# A representative cosmological set of base units (z ≈ 0 box), enough to make the
# unit factors finite and physical for the conservation checks.
units(T) = CK.chemistry_units(T; urho = 1.0e-28, uxyz = 3.086e24, utim = 3.156e15, aye = 1.0)

# Build a single-cell grid in the species-of-arrays layout.
function one_cell(T; nsp, dens, ge, ionized)
    names = CK.species_names(nsp)
    arr(v) = (a = zeros(T, 1); a[1] = v; a)
    # Neutral-ish or ionized primordial mix (mass densities, He = 0.24 by mass).
    X, Y = T(0.76), T(0.24)
    f = ionized ? T(0.999) : T(1e-4)
    nt = (; de   = arr(X * f * dens),
            HI   = arr(X * (1 - f) * dens),
            HII  = arr(X * f * dens),
            HeI  = arr(Y * (1 - f) * dens),
            HeII = arr(Y * f * dens),
            HeIII= arr(T(1e-10) * dens))
    if nsp ≥ 9
        nt = merge(nt, (; HM = arr(T(1e-12)*dens), H2I = arr(T(1e-6)*dens), H2II = arr(T(1e-14)*dens)))
    end
    if nsp ≥ 12
        nt = merge(nt, (; DI = arr(T(1e-5)*dens), DII = arr(T(1e-9)*dens), HDI = arr(T(1e-10)*dens)))
    end
    return CK.Species(nt), arr(ge), arr(dens)
end

@testset "ChemistryKernels" begin
    T = Float64
    u = units(T)
    rt = CK.build_rate_tables(T; nratec = 400, units = u, casebrates = false)

    @testset "rate tables finite & positive" begin
        @test all(isfinite, rt.k1) && all(rt.k1 .> 0)
        @test all(isfinite, rt.k2) && all(rt.k2 .> 0)
        @test all(isfinite, rt.reHII) && all(rt.reHII .> 0)
        @test rt.grid.nratec == 400
    end

    for nsp in (6, 9, 12)
        @testset "network nsp=$nsp conserves H, He, charge" begin
            sp, ge, dens = one_cell(T; nsp = nsp, dens = 1.0, ge = 1.0e3, ionized = true)
            # totals before
            Hbefore  = sp.HI[1] + sp.HII[1] + (nsp ≥ 9 ? sp.HM[1] + sp.H2I[1] + sp.H2II[1] : 0.0)
            Hebefore = sp.HeI[1] + sp.HeII[1] + sp.HeIII[1]
            CK.solve_rate_cool!(ge, dens, sp, rt, u;
                nspecies = nsp, gamma = 5/3, temperature_units = 1.0e4,
                dt = 1.0e-3, idim = 1, i1 = 1, i2 = 1, itmax = 5000)
            Hafter  = sp.HI[1] + sp.HII[1] + (nsp ≥ 9 ? sp.HM[1] + sp.H2I[1] + sp.H2II[1] : 0.0)
            Heafter = sp.HeI[1] + sp.HeII[1] + sp.HeIII[1]
            @test isapprox(Hafter, Hbefore; rtol = 1e-6)
            @test isapprox(Heafter, Hebefore; rtol = 1e-6)
            # charge conservation closure
            de_q = sp.HII[1] + 0.25*sp.HeII[1] + 0.5*sp.HeIII[1] +
                   (nsp ≥ 9 ? 0.5*sp.H2II[1] - sp.HM[1] : 0.0)
            @test isapprox(sp.de[1], max(de_q, 1e-20); rtol = 1e-6)
            @test all(>=(0), sp.HI) && all(>=(0), sp.HII)
            @test isfinite(ge[1]) && ge[1] > 0
            if nsp == 12
                @test isfinite(sp.DI[1]) && sp.DI[1] > 0
                @test isfinite(sp.DII[1]) && sp.DII[1] > 0
                @test isfinite(sp.HDI[1]) && sp.HDI[1] > 0
            end
        end
    end

    @testset "Enzo-parity mode (conserve=false) runs & stays finite" begin
        # bit-exact-Enzo mode skips make_consistent; species evolve by the raw
        # semi-implicit BDF + charge-conservation electrons, as solve_rate_cool.F.
        sp, ge, dens = one_cell(T; nsp = 9, dens = 1.0, ge = 1.0e3, ionized = false)
        CK.solve_rate_cool!(ge, dens, sp, rt, u;
            nspecies = 9, gamma = 5/3, temperature_units = 1.0e4,
            dt = 1.0e-3, idim = 1, i1 = 1, i2 = 1, itmax = 5000, conserve = false)
        @test all(isfinite, sp.HI) && all(isfinite, sp.H2I)
        @test sp.HI[1] > 0 && ge[1] > 0
    end

    @testset "hot gas cools (edot < 0)" begin
        sp, ge, dens = one_cell(T; nsp = 6, dens = 1.0, ge = 1.0e4, ionized = true)
        edot = zeros(T, 1)
        CK.cool1d_multi!(edot, ge, dens, sp, rt, u;
            nspecies = 6, gamma = 5/3, temperature_units = 1.0e4, idim = 1, i1 = 1, i2 = 1)
        @test edot[1] < 0          # collisionally ionized primordial gas radiates
    end
end
