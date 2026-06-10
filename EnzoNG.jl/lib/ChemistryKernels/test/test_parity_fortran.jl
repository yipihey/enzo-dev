# ── Fortran golden-fixture parity (chemistry) — RUNNABLE ──────────────────────
# Ground truth is Enzo's ACTUAL solve_rate_cool.F, compiled standalone into
# libenzochem.so by parity/build_enzo_chem.sh (gfortran; no full Enzo build).
# This test builds it (if gfortran + the src/enzo tree are present) and asserts
# ChemistryKernels.solve_rate_cool!(conserve=false) tracks Enzo. It self-skips
# cleanly otherwise, so it is safe to include from runtests.jl.
#
# Measured agreement (one-zone, dt-evolved): neutral fractions to ~1e-3; the
# ionized fraction in the thermally-unstable cooling-peak regime is sensitive
# (~15%); internal energy to a few %. Tolerances below reflect that.

using Test
using ChemistryKernels
const CK = ChemistryKernels
const mh = 1.67262171e-24; const kb = 1.3806504e-16

const PARITY_DIR = normpath(joinpath(@__DIR__, "..", "parity"))
const SRC_ENZO   = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "src", "enzo"))
const LIBPATH    = joinpath(PARITY_DIR, "build", "libenzochem.so")

function _build_enzochem()
    isfile(LIBPATH) && return true
    (Sys.which("gfortran") !== nothing && isdir(SRC_ENZO)) || return false
    try
        run(`bash $(joinpath(PARITY_DIR, "build_enzo_chem.sh"))`)
        return isfile(LIBPATH)
    catch err
        @info "enzochem build failed; skipping Fortran parity" err
        return false
    end
end

function enzo_step(s; ispecies, nratec, d, dt, units...)
    R(x) = Ref{Float64}(x); ierr = Ref{Int64}(0)
    refs = Dict(k => R(getfield(s, k)) for k in keys(s))
    u = (; units...)
    ccall(("enzo_chem_onezone_", LIBPATH), Cvoid,
        (Ref{Int64},Ref{Int64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
         Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
         Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
         Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
         Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
         Ref{Int64},Ref{Int64}),
        Ref{Int64}(ispecies), Ref{Int64}(nratec), R(d), refs[:ge], refs[:de], refs[:HI],
        refs[:HII], refs[:HeI], refs[:HeII], refs[:HeIII], refs[:HM], refs[:H2I],
        refs[:H2II], refs[:DI], refs[:DII], refs[:HDI], R(dt), R(u.dx), R(u.aye),
        R(0.0), R(u.temstart), R(u.temend), R(u.gamma), R(u.fh), R(u.dtoh),
        R(u.utem), R(u.uxyz), R(u.uaye), R(u.urho), R(u.utim), Ref{Int64}(0), ierr)
    return (; (k => refs[k][] for k in keys(s))...)
end

@testset "ChemistryKernels ⟷ Enzo solve_rate_cool parity" begin
    if !_build_enzochem()
        @test_skip "gfortran or src/enzo unavailable — Fortran parity skipped"
    else
        urho=1.0e-24; uxyz=3.086e21; utim=3.156e13; aye=1.0; uaye=1.0
        utem=(mh/kb)*(uxyz/utim)^2; gamma=5/3; fh=0.76; dtoh=3.4e-5
        temstart=1.0; temend=1.0e9; nratec=400; d=1.0
        un = (; dx=1.0, aye, temstart, temend, gamma, fh, dtoh, utem, uxyz, uaye, urho, utim)
        u  = CK.chemistry_units(Float64; urho, uxyz, utim, aye, uaye)
        rt = CK.build_rate_tables(Float64; nratec, temstart, temend, units=u, casebrates=false)

        for (tgas0, atol_HI, rtol_ge) in ((1.0e4, 2e-3, 2e-2), (3.0e4, 4e-3, 5e-2))
            ge0 = tgas0/((gamma-1)*1.22*utem)
            ic = (; ge=ge0, de=1e-4, HI=0.76, HII=1e-4, HeI=0.24, HeII=1e-6,
                   HeIII=1e-12, HM=1e-11, H2I=1e-6, H2II=1e-13, DI=dtoh*0.76,
                   DII=1e-10, HDI=1e-10)
            es = ic
            names = CK.species_names(6)
            nt = NamedTuple(n => [Float64(getfield(ic,n))] for n in names)
            sp = CK.Species(nt); ge=[ic.ge]; dens=[d]; dt=0.02
            for _ in 1:10
                es = enzo_step(es; ispecies=1, nratec, d, dt, un...)
                CK.solve_rate_cool!(ge, dens, sp, rt, u; nspecies=6, gamma,
                    temperature_units=utem, dt, idim=1, i1=1, i2=1, conserve=false, itmax=50000)
            end
            @testset "tgas0=$tgas0 K" begin
                @test isapprox(getproperty(sp,:HI)[1], es.HI; atol=atol_HI)   # neutral fraction
                @test isapprox(ge[1], es.ge; rtol=rtol_ge)                     # internal energy
                @test getproperty(sp,:HII)[1] > 0
            end
        end
    end
end
