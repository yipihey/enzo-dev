# Krome.jl — parser / expression-compiler / executor-bridge suite (CPU, headless).
# The network-file tests gate on the KROME submodule being checked out; the
# expression-compiler tests are self-contained.

using Test
using Krome
using ChemistryKernels
const CK = ChemistryKernels

@testset "Krome.jl" begin

    @testset "Fortran rate-expression translation" begin
        # constant rate, exercises the d-literal translation
        f1 = Krome.compile_rate("6.0d-10")
        @test f1 !== nothing
        @test f1(1234.0) ≈ 6.0e-10

        # Arrhenius in Tgas
        f2 = Krome.compile_rate("1.0d-9*exp(-4.1d1/Tgas)")
        @test f2(100.0) ≈ 1.0e-9 * exp(-41.0/100.0)

        # uses Te / invTe / ** → ^
        f3 = Krome.compile_rate("3.92d-13*invTe**0.6353d0")
        Te = 1.0e4 * 8.617343e-5
        @test f3(1.0e4) ≈ 3.92e-13 * (1/Te)^0.6353

        # density-dependent rate → unsupported → nothing (skipped cleanly)
        @test Krome.compile_rate("2.0d-10*n(idx_H)") === nothing
    end

    submods = Krome.list_krome_networks()
    if isempty(submods)
        @info "KROME submodule not checked out — skipping network-file tests. " *
              "Run: git submodule update --init EnzoNG.jl/lib/Krome/extern/krome"
    else
        @testset "parse react_primordial" begin
            net = Krome.parse_krome(Krome.krome_network_path("react_primordial"))
            @test "H" in net.species && "H+" in net.species && "E" in net.species
            @test "HE" in net.species
            @test length(net.reactions) ≥ 20
            # every kept reaction has a callable rate and valid indices
            for r in net.reactions
                @test all(s -> 1 ≤ s ≤ length(net.species), r.reactants)
                @test all(s -> 1 ≤ s ≤ length(net.species), r.products)
                @test isfinite(r.rate(1.0e3))
            end
        end

        @testset "lower to GenericNetwork and integrate on CPU" begin
            net = Krome.parse_krome(Krome.krome_network_path("react_primordial"))
            gen = Krome.to_generic(Float64, net; nratec = 200)
            @test gen.nspec == length(net.species)
            @test gen.nreac == length(net.reactions)

            be = CK.backend(:cpu)
            ncells = 4
            n = CK.device_zeros(be, Float64, (gen.nspec, ncells))
            # seed a warm, mostly-neutral primordial mix (cm^-3)
            si(name) = net.index[name]
            for c in 1:ncells
                n[si("H"), c]  = 1.0
                n[si("HE"), c] = 0.08
                n[si("E"), c]  = 1.0e-4
                n[si("H+"), c] = 1.0e-4
            end
            ntot0 = sum(n)
            Tgas = fill(8000.0, ncells)
            C, D = CK.allocate_workspace(be, Float64, gen.nspec, ncells)
            CK.generic_step!(n, Tgas, gen, C, D; dt = 1.0e7, itmax = 2000)
            @test all(isfinite, n)
            @test all(>=(0), n)
            # mass-action with charge-neutral seed: total nuclei stay O(1), no blow-up
            @test sum(n) < 5 * ntot0
        end
    end
end
