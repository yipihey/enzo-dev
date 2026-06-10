# Enzo-Fortran (ground truth) ⟷ ChemistryKernels (KA port): agreement over a
# dynamic ionization trajectory, plus a throughput benchmark.
#
#   FC=gfortran bash build_enzo_chem.sh           # build libenzochem.so first
#   julia --project=../test compare_enzo.jl
using Printf
using ChemistryKernels
const CK = ChemistryKernels
const LIB = joinpath(@__DIR__, "build", "libenzochem.so")
const mh = 1.67262171e-24; const kb = 1.3806504e-16

function enzo_onezone!(s; ispecies, nratec, d, dt, dx, aye, redshift, temstart,
                       temend, gamma, fh, dtoh, utem, uxyz, uaye, urho, utim)
    R(x)=Ref{Float64}(x); ierr=Ref{Int64}(0)
    ge=R(s.ge); de=R(s.de); HI=R(s.HI); HII=R(s.HII); HeI=R(s.HeI); HeII=R(s.HeII)
    HeIII=R(s.HeIII); HM=R(s.HM); H2I=R(s.H2I); H2II=R(s.H2II); DI=R(s.DI)
    DII=R(s.DII); HDI=R(s.HDI)
    ccall(("enzo_chem_onezone_", LIB), Cvoid,
        (Ref{Int64},Ref{Int64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
         Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
         Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
         Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
         Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
         Ref{Int64},Ref{Int64}),
        Ref{Int64}(ispecies), Ref{Int64}(nratec), R(d), ge, de, HI, HII, HeI, HeII,
        HeIII, HM, H2I, H2II, DI, DII, HDI, R(dt), R(dx), R(aye), R(redshift),
        R(temstart), R(temend), R(gamma), R(fh), R(dtoh), R(utem), R(uxyz),
        R(uaye), R(urho), R(utim), Ref{Int64}(0), ierr)
    return (; ge=ge[], de=de[], HI=HI[], HII=HII[], HeI=HeI[], HeII=HeII[],
            HeIII=HeIII[], HM=HM[], H2I=H2I[], H2II=H2II[], DI=DI[], DII=DII[], HDI=HDI[])
end

# shared units & a HOT, mostly-neutral IC so the gas actively ionizes
urho=1.0e-24; uxyz=3.086e21; utim=3.156e13; aye=1.0; uaye=1.0
utem=(mh/kb)*(uxyz/utim)^2
gamma=5/3; fh=0.76; dtoh=3.4e-5; temstart=1.0; temend=1.0e9; nratec=400
d=1.0; tgas0=3.0e4
ge0=tgas0/((gamma-1)*1.22*utem)
ic = (; ge=ge0, de=1e-4, HI=0.76, HII=1e-4, HeI=0.24, HeII=1e-6, HeIII=1e-12,
       HM=1e-11, H2I=1e-6, H2II=1e-13, DI=dtoh*0.76, DII=1e-10, HDI=1e-10)

isfeasible = isfile(LIB)
isfeasible || error("build libenzochem.so first: FC=gfortran bash build_enzo_chem.sh")

println("# Enzo ⟷ ChemistryKernels parity  (utem=$(round(utem)) K, tgas0=$tgas0 K)\n")
for ispecies in (1,2,3)
    nsp=(6,9,12)[ispecies]
    names=CK.species_names(nsp)
    u=CK.chemistry_units(Float64; urho, uxyz, utim, aye, uaye)
    rt=CK.build_rate_tables(Float64; nratec, temstart, temend, units=u, casebrates=false)
    # evolve both independently from the same IC over N macro-steps
    es=ic; sp_vals=Dict(pairs(ic)...)
    nt=NamedTuple(n=>[Float64(getfield(ic,n))] for n in names); sp=CK.Species(nt)
    ge=[ic.ge]; dens=[d]; dt=0.02; N=10
    maxre=Dict(f=>0.0 for f in (:HI,:HII,:HeII,:HeIII,:H2I,:de,:ge))
    for step in 1:N
        es=enzo_onezone!(es; ispecies, nratec, d, dt, dx=1.0, aye, redshift=0.0,
            temstart, temend, gamma, fh, dtoh, utem, uxyz, uaye, urho, utim)
        CK.solve_rate_cool!(ge, dens, sp, rt, u; nspecies=nsp, gamma,
            temperature_units=utem, dt, idim=1, i1=1, i2=1, conserve=false, itmax=50000)
        for f in keys(maxre)
            (f in names || f===:ge) || continue
            ev=getfield(es,f); cv=(f===:ge) ? ge[1] : getproperty(sp,f)[1]
            maxre[f]=max(maxre[f], abs(cv-ev)/max(abs(ev),1e-30))
        end
    end
    @printf "ispecies=%d (%2d-sp): HI=%.6e/%.6e  HII=%.3e/%.3e (Enzo/CK)\n" ispecies nsp es.HI getproperty(sp,:HI)[1] es.HII getproperty(sp,:HII)[1]
    print("   max rel.err over trajectory: ")
    for f in (:HI,:HII,:HeII,:HeIII,:de,:ge); print(f,"=",@sprintf("%.1e ",maxre[f])); end
    nsp≥9 && print("H2I=",@sprintf("%.1e",maxre[:H2I]))
    println()
end

# ── throughput: many independent cells, one macro-step (CK on CPU) ───────────
println("\n# Throughput (ChemistryKernels, CPU, 9-species)")
u=CK.chemistry_units(Float64; urho, uxyz, utim, aye, uaye)
rt=CK.build_rate_tables(Float64; nratec, temstart, temend, units=u, casebrates=false)
for ncell in (10^3, 10^4, 10^5)
    names=CK.species_names(9)
    nt=NamedTuple(n=>fill(Float64(getfield(ic,n)),ncell) for n in names)
    sp=CK.Species(nt); ge=fill(ic.ge,ncell); dens=fill(d,ncell)
    CK.solve_rate_cool!(ge,dens,sp,rt,u; nspecies=9, gamma, temperature_units=utem,
        dt=0.02, idim=ncell, i1=1, i2=ncell, conserve=false, itmax=50000) # warm up
    t=@elapsed CK.solve_rate_cool!(ge,dens,sp,rt,u; nspecies=9, gamma,
        temperature_units=utem, dt=0.02, idim=ncell, i1=1, i2=ncell, conserve=false, itmax=50000)
    @printf "  %7d cells: %8.2f ms  (%6.2f µs/cell)\n" ncell t*1e3 t/ncell*1e6
end
