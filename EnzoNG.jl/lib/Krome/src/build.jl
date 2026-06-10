# ── KromeNetwork → ChemistryKernels.GenericNetwork ────────────────────────────
# Lowers a parsed KROME network onto the backend-agnostic mass-action executor in
# ChemistryKernels, tabulating each reaction's compiled `Tgas→k` closure on the
# shared log(T) grid. The result runs on CPU and GPU through KernelAbstractions.

"""
    to_generic(T=Float64, net::KromeNetwork; nratec=400, temstart=1.0, temend=1e9)
        -> ChemistryKernels.GenericNetwork

Tabulate `net`'s reaction rates and assemble a `GenericNetwork` of element type
`T`. Reaction `r`'s coefficient is `net.reactions[r].rate(Tgas)` (cgs). Use
`ChemistryKernels.generic_step!` to integrate; species are number densities (cm⁻³)
and `dt` is in seconds (KROME's cgs convention).
"""
function to_generic(::Type{T}, net::KromeNetwork; nratec::Integer = 400,
                    temstart = 1.0, temend = 1.0e9) where {T}
    reactions = [(r.reactants, r.products) for r in net.reactions]
    # Pre-tabulate the compiled rate closures (each already dispatches through
    # invokelatest, so direct calls here are world-age-safe), then hand the
    # finished table to the executor's precomputed-table path.
    nreac = length(net.reactions)
    l0 = log(T(temstart)); dl = (log(T(temend)) - l0) / (nratec - 1)
    tbl = Matrix{T}(undef, nratec, nreac)
    for r in 1:nreac, i in 1:nratec
        Tgas = exp(l0 + (i - 1) * dl)
        tbl[i, r] = T(net.reactions[r].rate(Tgas))
    end
    return ChemistryKernels.generic_network(T;
        nspec = length(net.species), reactions = reactions, rate_table = tbl,
        nratec = nratec, temstart = temstart, temend = temend)
end
to_generic(net::KromeNetwork; kw...) = to_generic(Float64, net; kw...)
