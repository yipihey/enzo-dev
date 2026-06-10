# ── KromeNetwork → executors ──────────────────────────────────────────────────
# Two lowering paths:
#   • to_generic  — tabulate the TEMPERATURE-ONLY reactions onto
#     ChemistryKernels.GenericNetwork (CPU+GPU, KernelAbstractions). Density-
#     dependent reactions can't be tabulated by T alone and are dropped here.
#   • direct_step! — a CPU mass-action integrator that evaluates the compiled
#     rate closures directly, so it runs the FULL network including the
#     density-dependent (ntot/Hnuclei-bridged) rates of the CO/metal networks.

"""
    to_generic(T=Float64, net::KromeNetwork; nratec=400, temstart=1.0, temend=1e9)
        -> (GenericNetwork, ndropped)

Tabulate the temperature-only reactions of `net` onto a GPU-capable
`ChemistryKernels.GenericNetwork`. Returns the network and the number of
density-dependent reactions dropped (run those with [`direct_step!`](@ref)).
"""
function to_generic(::Type{T}, net::KromeNetwork; nratec::Integer = 400,
                    temstart = 1.0, temend = 1.0e9) where {T}
    keep = [r for r in net.reactions if !r.density_dep]
    ndropped = length(net.reactions) - length(keep)
    reactions = [(r.reactants, r.products) for r in keep]
    nreac = length(keep)
    l0 = log(T(temstart)); dl = (log(T(temend)) - l0) / (nratec - 1)
    tbl = Matrix{T}(undef, nratec, nreac)
    for r in 1:nreac, i in 1:nratec
        Tgas = exp(l0 + (i - 1) * dl)
        tbl[i, r] = T(keep[r].rate(Tgas, 1.0))     # ntot unused for T-only rates
    end
    gen = ChemistryKernels.generic_network(T; nspec = length(net.species),
        reactions = reactions, rate_table = tbl, nratec = nratec,
        temstart = temstart, temend = temend)
    return gen, ndropped
end
to_generic(net::KromeNetwork; kw...) = to_generic(Float64, net; kw...)

"""
    direct_step!(n, Tgas, net::KromeNetwork; dt, itmax=10000, ntot=nothing) -> n

CPU mass-action integrator over the FULL `net` (including density-dependent
rates), advancing the `nspec×ncells` number-density matrix `n` (cm⁻³) by `dt`
seconds at per-cell temperatures `Tgas`. Per-cell total density `ntot` (for the
`Hnuclei`/`ntot`-bridged rates) defaults to the column sum of `n`. Uses the same
Anninos+97 semi-implicit backward-difference scheme as the primordial solver.
"""
function direct_step!(n::AbstractMatrix{T}, Tgas::AbstractVector, net::KromeNetwork;
                      dt::Real, itmax::Integer = 10000, ntot = nothing) where {T}
    nspec, ncells = size(n)
    rxn = net.reactions; nreac = length(rxn)
    C = zeros(T, nspec); D = zeros(T, nspec)
    for c in 1:ncells
        Tc = T(Tgas[c]); ttot = zero(T)
        for _ in 1:itmax
            nt = ntot === nothing ? sum(@view n[:, c]) : T(ntot[c])
            fill!(C, zero(T)); fill!(D, zero(T))
            for r in 1:nreac
                rr = rxn[r]
                k = T(rr.rate(Tc, nt))
                R = k
                for s in rr.reactants; R *= n[s, c]; end
                for p in rr.products; C[p] += R; end
                for s in rr.reactants; D[s] += R / max(n[s, c], T(1e-30)); end
            end
            dtsub = dt - ttot
            for s in 1:nspec
                D[s] > 1e-30 && (dtsub = min(dtsub, T(0.1) / D[s]))
            end
            dtsub = max(dtsub, T(1e-6) * dt)
            for s in 1:nspec
                n[s, c] = max((C[s] * dtsub + n[s, c]) / (one(T) + D[s] * dtsub), T(1e-30))
            end
            ttot += dtsub
            ttot ≥ dt * (one(T) - T(1e-6)) && break
        end
    end
    return n
end
