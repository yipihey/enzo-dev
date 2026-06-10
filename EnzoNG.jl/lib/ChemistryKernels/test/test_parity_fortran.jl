# ── Fortran golden-fixture parity scaffold (chemistry) ────────────────────────
# The certification layer for ChemistryKernels, mirroring PPMKernels/test/harness:
#
#   Layer A  port faithfulness : KA solve_rate_cool! on CPU/Float64  vs  Enzo's
#            Fortran solve_rate_cool (via EnzoLib's session bridge), bit-tight.
#   Layer C  f32 accuracy floor : Metal/Float32  vs  the f64 Fortran reference.
#
# Enzo's native scheme does NOT renormalize species (charge-conservation electrons
# only), so parity runs the KA solver with `conserve=false`. The Fortran reference
# is `enzomodules_session_solve_cooling` (Grid_SolveRadiativeCooling → multi_cool
# → solve_rate_cool.F) over a one-zone multispecies grid.
#
# This file is GATED: it self-skips unless the live-Enzo bridge is built and
# loadable (EnzoLib.grid_available()). It is NOT part of the headless runtests.jl;
# run it explicitly once the bridge exists:
#     <julia> --project=test test/test_parity_fortran.jl
# See docs/adr/0006 for the build + the remaining field-accessor hookup.

using Test
using ChemistryKernels
const CK = ChemistryKernels

const _ENZO = try
    @eval using EnzoLib
    EnzoLib.grid_available()
catch err
    @info "EnzoLib bridge unavailable — chemistry parity layer skips" err
    false
end

@testset "ChemistryKernels ⟷ Enzo solve_rate_cool parity" begin
    if !_ENZO
        @test_skip "Enzo bridge not built (EnzoLib.grid_available() == false)"
    else
        # Intended comparison (one-zone, MultiSpecies networks 1/2/3):
        #   1. stage identical (de,HI,…,HDI, ge, ρ) into an Enzo one-zone grid and
        #      a matching CK.Species on the CPU backend;
        #   2. step Enzo:  EnzoLib.session_solve_cooling(h, level)  (Fortran);
        #   3. step CK:     CK.solve_rate_cool!(…; nspecies, dt, conserve=false);
        #   4. read Enzo's post-step fields back and diff against the CK fields
        #      under Layer-A tolerance (rtol≈1e-12), and the Metal/f32 run under
        #      Layer-C (rtol≈1e-5).
        #
        # The remaining hookup is a one-zone field accessor on the bridge
        # (stage/extract the species + energy arrays); the rate/cooling tables are
        # built from the SAME calc_rates inputs Enzo used, so Layer A is bit-tight
        # by construction once the fields are wired. Tracked in ADR-0006.
        @test_broken false   # placeholder until the one-zone field accessor lands
    end
end
