# KROME ⟷ Enzo chemistry: agreement (rate-level + collisional-ionization
# equilibrium) and performance. KROME's react_primordial and Enzo's
# calc_rates.F are both Abel+97-based, so the rate sets should agree closely; the
# generic mass-action executor (KROME path) and the specialised solve_rate_cool!
# (Enzo port) should reach the same ionization equilibrium.
#
#   julia --project=../test compare_krome_enzo.jl
using Printf
using Krome
using ChemistryKernels
const CK = ChemistryKernels
const mh = 1.67262171e-24; const kb = 1.3806504e-16

# --- build both networks ----------------------------------------------------
net = Krome.parse_krome(Krome.krome_network_path("react_primordial"))
gen, ndrop = Krome.to_generic(Float64, net; nratec = 400)
@printf "KROME react_primordial: %d species, %d reactions (%d tabulated, %d density-dep dropped, %d skipped)\n" length(net.species) length(net.reactions) gen.nreac ndrop net.skipped

# Enzo rate tables (cgs via ×kunit) for the same conditions
urho=1.0e-24; uxyz=3.086e21; utim=3.156e13; aye=1.0; uaye=1.0
u  = CK.chemistry_units(Float64; urho, uxyz, utim, aye, uaye)
rt = CK.build_rate_tables(Float64; nratec=400, temstart=1.0, temend=1.0e9, units=u, casebrates=false)
kunit = u.kunit
enzo_k(tbl, T) = (g=rt.grid; lt=log(clamp(T,1.0,1e9)); i=clamp(Int(floor((lt-g.logtem0)/g.dlogtem))+1,1,g.nratec-1); td=(lt-(g.logtem0+(i-1)*g.dlogtem))/g.dlogtem; (tbl[i]+td*(tbl[i+1]-tbl[i]))*kunit)

# find KROME reactions by species signature
sset(x) = Set(net.index[s] for s in x)
function find_rxn(reac, prod)
    for r in net.reactions
        Set(r.reactants)==sset(reac) && Set(r.products)==sset(prod) && return r
    end
    nothing
end
kHIion = find_rxn(["H","E"], ["H+","E","E"])     # HI + e -> HII + 2e   (Enzo k1)
# recombination is piecewise across T branches — SUM all matching reactions
sumrate(reac, T) = sum(r.rate(T,1.0) for r in net.reactions if Set(r.reactants)==sset(reac); init=0.0)

println("\n# Rate agreement  (cgs cm^3/s):  KROME vs Enzo(calc_rates)")
@printf "%8s | %-22s | %-22s\n" "T [K]" "HI ionization (k1)" "HII recomb (k2)"
for T in (1.0e4, 3.0e4, 1.0e5, 3.0e5, 1.0e6)
    kK1 = kHIion.rate(T, 1.0); kE1 = enzo_k(rt.k1, T)
    kK2 = sumrate(["H+","E"], T); kE2 = enzo_k(rt.k2, T)
    @printf "%8.0e | K %.3e E %.3e (%4.1f%%) | K %.3e E %.3e (%4.1f%%)\n" T kK1 kE1 100*abs(kK1-kE1)/kE1 kK2 kE2 100*abs(kK2-kE2)/kE2
end

# --- collisional ionization equilibrium (isothermal) ------------------------
# Evolve the KROME network at fixed T to equilibrium; compare HII/H_tot to the
# Enzo two-level estimate k1/(k1+k2) (the dominant H ionization balance).
println("\n# Collisional ionization equilibrium:  HII / H_total")
@printf "%8s | %-12s | %-12s | %-8s\n" "T [K]" "KROME (gen)" "Enzo k1/(k1+k2)" "rel.diff"
be = CK.backend(:cpu)
iH=net.index["H"]; iHp=net.index["H+"]; iE=net.index["E"]
iHe=net.index["HE"]; iHep=net.index["HE+"]; iHepp=net.index["HE++"]
for T in (1.5e4, 2.0e4, 3.0e4, 5.0e4, 1.0e5)
    n = zeros(Float64, gen.nspec, 1)
    n[iH,1]=1.0; n[iHe,1]=0.0; n[iE,1]=1e-6; n[iHp,1]=1e-6
    Tg=fill(T,1); C,D=CK.allocate_workspace(be,Float64,gen.nspec,1)
    for _ in 1:400  # march to equilibrium (longer at low T where ionization is slow)
        CK.generic_step!(n, Tg, gen, C, D; dt=3.0e13, itmax=20000)
    end
    xK = n[iHp,1]/(n[iH,1]+n[iHp,1])
    k1=enzo_k(rt.k1,T)/kunit; k2=enzo_k(rt.k2,T)/kunit  # ratio is unit-independent
    xE = k1/(k1+k2)
    @printf "%8.0e | %-12.4f | %-12.4f | %6.1f%%\n" T xK xE 100*abs(xK-xE)/xE
end

# --- performance ------------------------------------------------------------
println("\n# Throughput (CPU, per cell, one macro-step)")
names9 = CK.species_names(9)
ic9 = (de=1e-4,HI=0.76,HII=1e-4,HeI=0.24,HeII=1e-6,HeIII=1e-12,HM=1e-11,H2I=1e-6,H2II=1e-13)
utem=(mh/kb)*(uxyz/utim)^2; ge0=3e4/((5/3-1)*1.22*utem)
for ncell in (10^4, 10^5)
    # Enzo-port specialised solver
    nt=NamedTuple(s=>fill(Float64(getfield(ic9,s)),ncell) for s in names9)
    sp=CK.Species(nt); ge=fill(ge0,ncell); dens=fill(1.0,ncell)
    CK.solve_rate_cool!(ge,dens,sp,rt,u; nspecies=9,gamma=5/3,temperature_units=utem,dt=0.02,idim=ncell,i1=1,i2=ncell,conserve=false)
    tCK=@elapsed CK.solve_rate_cool!(ge,dens,sp,rt,u; nspecies=9,gamma=5/3,temperature_units=utem,dt=0.02,idim=ncell,i1=1,i2=ncell,conserve=false)
    # KROME generic network
    ng=zeros(Float64,gen.nspec,ncell); ng[iH,:].=0.76; ng[iHp,:].=1e-4; ng[iE,:].=1e-4; ng[iHe,:].=0.08
    Tg=fill(3.0e4,ncell); C,D=CK.allocate_workspace(be,Float64,gen.nspec,ncell)
    CK.generic_step!(ng,Tg,gen,C,D; dt=6.3e11, itmax=2000)
    tK=@elapsed CK.generic_step!(ng,Tg,gen,C,D; dt=6.3e11, itmax=2000)
    @printf "  %7d cells:  Enzo-port solve_rate_cool! %.2f µs/cell   KROME generic_step! %.2f µs/cell\n" ncell tCK/ncell*1e6 tK/ncell*1e6
end
