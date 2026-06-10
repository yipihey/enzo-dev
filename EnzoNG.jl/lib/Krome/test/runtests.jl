# Krome.jl — parser / expression-compiler / executor-bridge suite (CPU, headless).
# Network-file tests gate on the KROME submodule; expression tests are standalone.

using Test
using Krome
using ChemistryKernels
const CK = ChemistryKernels

# compile_rate now returns ((Tgas,ntot)->k, density_dep) or nothing
rate1(s; kw...) = (c = Krome.compile_rate(s; kw...); c === nothing ? nothing : c[1])

@testset "Krome.jl" begin

    @testset "Fortran rate-expression translation" begin
        f1 = rate1("6.0d-10")
        @test f1 !== nothing && f1(1234.0, 1.0) ≈ 6.0e-10
        f2 = rate1("1.0d-9*exp(-4.1d1/Tgas)")
        @test f2(100.0, 1.0) ≈ 1.0e-9 * exp(-41.0/100.0)
        f3 = rate1("3.92d-13*invTe**0.6353d0")
        Te = 1.0e4 * 8.617343e-5
        @test f3(1.0e4, 1.0) ≈ 3.92e-13 * (1/Te)^0.6353
        # density-dependent rate now compiles via @var bridging + ntot
        c = Krome.compile_rate("kl11^(1.-a11)*kh11^a11";
            vars = [:logT=>"log10(Tgas)", :invT=>"1d0/Tgas",
                    :kl11=>"1d1**(-27.029d0+3.801d0*logT-29487d0*invT)",
                    :kh11=>"1d1**(-2.729d0-1.75d0*logT-23474d0*invT)",
                    :ncr11=>"1d1**(5.0792d0*(1d0-1.23d-5*(Tgas-2d3)))",
                    :a11=>"1.d0/(1.d0+(Hnuclei/(ncr11+1d-40)))"])
        @test c !== nothing
        f, dd = c
        @test dd                                    # flagged density-dependent
        @test f(2000.0, 1e2) != f(2000.0, 1e8)      # varies with ntot
        # genuinely unsupported (explicit species density) → skipped
        @test Krome.compile_rate("2.0d-10*n(idx_C)") === nothing
    end

    submods = Krome.list_krome_networks()
    if isempty(submods)
        @info "KROME submodule not checked out — skipping network-file tests."
    else
        @testset "parse react_primordial → GenericNetwork → integrate (CPU)" begin
            net = Krome.parse_krome(Krome.krome_network_path("react_primordial"))
            @test "H" in net.species && "H+" in net.species && "E" in net.species
            @test length(net.reactions) ≥ 20
            for r in net.reactions
                @test isfinite(r.rate(1.0e3, 1.0))
            end
            gen, ndrop = Krome.to_generic(Float64, net; nratec = 200)
            @test gen.nspec == length(net.species)
            @test gen.nreac == count(r -> !r.density_dep, net.reactions)
            be = CK.backend(:cpu); ncell = 4
            n = CK.device_zeros(be, Float64, (gen.nspec, ncell))
            si(x) = net.index[x]
            for c in 1:ncell
                n[si("H"),c]=1.0; n[si("HE"),c]=0.08; n[si("E"),c]=1e-4; n[si("H+"),c]=1e-4
            end
            ntot0 = sum(n); Tgas = fill(8000.0, ncell)
            Cw, Dw = CK.allocate_workspace(be, Float64, gen.nspec, ncell)
            CK.generic_step!(n, Tgas, gen, Cw, Dw; dt = 1.0e7, itmax = 2000)
            @test all(isfinite, n) && all(>=(0), n) && sum(n) < 5*ntot0
        end

        @testset "density-dependent CO/metal network compiles & runs (CPU direct)" begin
            net = Krome.parse_krome(Krome.krome_network_path("react_COthin"))
            ndd = count(r -> r.density_dep, net.reactions)
            @info "react_COthin" species=length(net.species) reactions=length(net.reactions) density_dependent=ndd skipped=net.skipped
            @test length(net.reactions) > 50
            @test ndd ≥ 1                                  # @var/ntot bridging worked
            # run the FULL network (incl. density-dependent) on CPU
            ncell = 2; n = zeros(Float64, length(net.species), ncell)
            seed(x,v) = haskey(net.index,x) && (n[net.index[x],:] .= v)
            seed("H",1e2); seed("H2",1e2); seed("HE",1e1); seed("E",1e-2)
            seed("C",1e-2); seed("O",1e-2); seed("H+",1e-2)
            Tgas = fill(50.0, ncell)
            Krome.direct_step!(n, Tgas, net; dt = 1.0e10, itmax = 2000)
            @test all(isfinite, n) && all(>=(0), n)
        end
    end
end
