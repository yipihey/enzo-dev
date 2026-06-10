# Grackle.jl — Cloudy-table loader + KA metal-cooling + libgrackle binding.
# Table tests gate on a Cloudy HDF5 file in the enzo-dev tree; the libgrackle
# binding test gates on a built library (GRACKLE_LIB / deps/build).

using Test
using Grackle
using ChemistryKernels
const CK = ChemistryKernels

# a Cloudy metals table shipped with enzo-dev (rank-3, sliced at z=0)
const TBL = let root = normpath(joinpath(@__DIR__, "..", "..", "..", "..")),
                cand = joinpath(root, "run", "Cooling", "CoolingTest_Cloudy", "solar_2008_3D_metals.h5")
    isfile(cand) ? cand : ""
end

@testset "Grackle.jl" begin

    if isempty(TBL)
        @info "No Cloudy HDF5 table found in tree — skipping table tests."
    else
        tbl = Grackle.load_cloudy_table(TBL; redshift = 0.0)

        @testset "load + grid sanity" begin
            @test tbl.rank == 2
            @test tbl.n1 > 1 && tbl.n2 > 1
            @test tbl.par1[1] < tbl.par1[end]      # log n_H ascending
            @test tbl.par2[1] < tbl.par2[end]      # log T ascending
        end

        @testset "KA interpolation = grid-node values (port fidelity)" begin
            # interpolating exactly at a grid node returns that node's value
            for (i1, i2) in ((5, 10), (20, 40), (30, 60))
                x1 = tbl.par1[i1]; x2 = tbl.par2[i2]
                node = tbl.cooling[(i1 - 1) * tbl.n2 + i2]
                @test CK._interp2d(x1, x2, tbl, tbl.cooling) ≈ node rtol=1e-12
            end
        end

        @testset "KA Cloudy metal cooling: physical shape" begin
            Ts = [1e4, 3e4, 1e5, 3e5, 1e6, 1e7]
            rhoH = fill(1.0, length(Ts)); Z = fill(1.0, length(Ts))
            em = zeros(length(Ts))
            CK.cloudy_metal_cooling!(em, rhoH, log.(Ts), Z, tbl;
                dom = 0.76, idim = length(Ts), i1 = 1, i2 = length(Ts), zscale = true)
            @test all(isfinite, em) && all(<(0), em)           # all cooling
            @test argmin(em) in (3, 4)                         # peaks at 1e5–3e5 K (metal lines)
            # cooling scales linearly with metallicity
            em2 = zeros(1)
            CK.cloudy_metal_cooling!(em2, [1.0], [log(1e5)], [2.0], tbl;
                dom = 0.76, idim = 1, i1 = 1, i2 = 1, zscale = true)
            @test em2[1] ≈ 2 * em[3] rtol=1e-12
        end

        @testset "integration: metals shorten the cooling time" begin
            mh = 1.67262171e-24; kb = 1.3806504e-16
            urho=1.0e-24; uxyz=3.086e21; utim=3.156e13
            utem=(mh/kb)*(uxyz/utim)^2
            u  = CK.chemistry_units(Float64; urho, uxyz, utim, aye=1.0)
            rt = CK.build_rate_tables(Float64; nratec=400, temstart=1.0, temend=1e9, units=u)
            ge0 = 1e5/((5/3-1)*0.6*utem)
            mk() = (nt = NamedTuple(s => [Float64(v)] for (s,v) in
                        (:de=>1e-1,:HI=>1e-3,:HII=>0.76,:HeI=>1e-4,:HeII=>0.24,:HeIII=>1e-4));
                    CK.Species(nt))
            # Z = 0 (primordial only)
            sp0 = mk(); ge0a=[ge0]; d=[1.0]
            CK.solve_rate_cool!(ge0a, d, sp0, rt, u; nspecies=6, gamma=5/3,
                temperature_units=utem, dt=0.01, idim=1, i1=1, i2=1, conserve=false)
            # Z = 1 with Cloudy metals
            sp1 = mk(); ge1a=[ge0]
            CK.solve_rate_cool!(ge1a, d, sp1, rt, u; nspecies=6, gamma=5/3,
                temperature_units=utem, dt=0.01, idim=1, i1=1, i2=1, conserve=false,
                cloudy=tbl, metallicity=1.0, fh=0.76)
            @test ge1a[1] < ge0a[1]      # metal cooling removes MORE energy
        end
    end

    @testset "libgrackle binding" begin
        if Grackle.grackle_available()
            v = Grackle.grackle_version()
            @info "libgrackle" version=v.version branch=v.branch
            @test !isempty(v.version)
        else
            @test_skip "libgrackle not built (set GRACKLE_LIB or run deps/build_grackle.sh)"
        end
    end
end
