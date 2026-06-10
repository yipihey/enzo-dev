# RadiationKernels — CPU correctness suite for the M1 moment solver.
#
# Headless f64-on-CPU checks of the closure limits, free-streaming transport, and
# the chemistry coupling. The bit-tight parity layer (vs an Enzo M1/ray-trace
# reference) slots in alongside, mirroring PPMKernels/test.

using Test
using RadiationKernels
const RK = RadiationKernels

@testset "RadiationKernels" begin

    @testset "M1 Eddington factor limits" begin
        @test RK.eddington_factor(0.0) ≈ 1/3       # isotropic
        @test RK.eddington_factor(1.0) ≈ 1.0       # free streaming
        @test 1/3 ≤ RK.eddington_factor(0.5) ≤ 1.0 # monotone in between
        @test RK.eddington_factor(-0.1) ≈ 1/3      # clamped
        @test RK.eddington_factor(1.5) ≈ 1.0       # clamped
    end

    @testset "free-streaming pulse advances at ~c and conserves photons" begin
        T = Float64
        idim, jdim, kdim = 128, 1, 1
        be = RK.backend(:cpu)
        c  = RK.SPEED_OF_LIGHT
        dx = 1.0
        dt = dx / (1 * c) * 0.9            # 1-D split CFL
        fld  = RK.RadField(be, T, (idim*jdim*kdim,))
        scr  = RK.RadField(be, T, (idim*jdim*kdim,))
        # a right-moving beam well inside the domain: f = |F|/(cN) = 1, fully
        # forward-peaked. Kept away from both edges so the (transmissive) boundary
        # faces stay at N≈0 and the conservative telescoping sum is clean.
        for i in 40:56
            fld.N[i]  = 1.0
            fld.Fx[i] = c * 1.0
        end
        N0 = sum(fld.N)
        x0 = sum((1:idim) .* fld.N) / N0  # photon centroid
        nsteps = 12
        for _ in 1:nsteps
            RK.m1_step!(fld, scr, c, dt, dx; idim=idim, jdim=jdim, kdim=kdim, ndim=1)
        end
        N1 = sum(fld.N)
        x1 = sum((1:idim) .* fld.N) / N1
        @test isapprox(N1, N0; rtol = 1e-3)           # photon number conserved
        @test x1 > x0 + 4                              # advected to the right (GLF diffuses)
        @test all(>=(0), fld.N)                        # positivity
    end

    @testset "photo rates accumulate over groups & feed chemistry" begin
        T = Float64
        be = RK.backend(:cpu)
        nd = 8
        g1 = RK.RadField(be, T, (nd,)); fill!(g1.N, 1.0e-3)
        g2 = RK.RadField(be, T, (nd,)); fill!(g2.N, 2.0e-3)
        xs = [RK.PhotoCrossSections(6.3e-18, 1.0e-11, 0.0, 0.0, 0.0, 0.0, 0.0),
              RK.PhotoCrossSections(1.2e-18, 3.0e-11, 0.0, 0.0, 0.0, 0.0, 0.0)]
        nHI  = fill(1.0, nd); nHeI = zeros(nd); nHeII = zeros(nd)
        k = (zeros(T, nd), zeros(T, nd), zeros(T, nd), zeros(T, nd), zeros(T, nd))
        RK.photo_rates!(k..., [g1, g2], xs, nHI, nHeI, nHeII; c = RK.SPEED_OF_LIGHT)
        expected = RK.SPEED_OF_LIGHT * (6.3e-18*1.0e-3 + 1.2e-18*2.0e-3)
        @test isapprox(k[1][1], expected; rtol = 1e-10)   # kphHI summed over groups
        @test k[5][1] > 0                                  # photo-heating positive
    end
end
