# ── Generic mass-action network executor ──────────────────────────────────────
# The hand-written `solve_rate_cool!` above is the certified primordial network.
# This is its *general* sibling: a precision-generic, backend-agnostic integrator
# for an ARBITRARY reaction network — the target the `Krome` package generates
# into, so any KROME `react_*` file becomes a CPU+GPU chemistry solver.
#
# A network is mass-action: reaction r with rate coefficient k_r(T) and up to 3
# reactants consumes them and produces up to 4 products. Per cell we advance the
# species with the same Anninos+97 semi-implicit backward-difference scheme the
# primordial solver uses, generalised:
#
#     Cₛ = Σ_{r: s a product}  (mult)·k_r·∏ n_reactant       (creation)
#     Dₛ = Σ_{r: s a reactant} (mult)·k_r·∏ n_other-reactant (destruction coef)
#     nₛ_new = (Cₛ·dt + nₛ) / (1 + Dₛ·dt)
#
# sub-cycled per cell. Species live in an `NS×ncells` column-major matrix (each
# cell's species contiguous); the per-cell creation/destruction accumulators are
# device scratch of the same shape. Temperature is supplied per cell (energy
# coupling / cooling is network-specific and layered on top, as in KROME).

export GenericNetwork, generic_network, generic_step!, allocate_workspace

"""
    GenericNetwork{IM,RM,T}

A tabulated mass-action reaction network ready for [`generic_step!`](@ref):

  * `nspec`, `nreac` — species / reaction counts.
  * `reac::IM` — `3×nreac` reactant species indices (column-major), `0` = empty slot.
  * `prod::IM` — `4×nreac` product species indices, `0` = empty slot.
  * `rate::RM` — `nratec×nreac` rate-coefficient table vs log(T), one column/reaction.
  * `grid::TempGrid{T}` — the shared log(T) tabulation grid.

Build with [`generic_network`](@ref) (host) and move to a device with
`adapt`/`to_device`-style helpers via the package's Adapt rules.
"""
struct GenericNetwork{IM,RM,T}
    nspec::Int
    nreac::Int
    reac::IM
    prod::IM
    rate::RM
    grid::TempGrid{T}
end

"""
    generic_network(T=Float64; nspec, reactions, rate_of, nratec=400,
                    temstart=1.0, temend=1e9) -> GenericNetwork{...}

Assemble a `GenericNetwork` from a list of `reactions`, each a
`(reactants::Vector{Int}, products::Vector{Int})` of 1-based species indices
(≤3 reactants, ≤4 products), and `rate_of(r, Tgas)::Real` giving reaction `r`'s
coefficient (already in the working code units) at temperature `Tgas`. The rates
are tabulated on a uniform log(T) grid, exactly like `calc_rates`.
"""
function generic_network(::Type{T} = Float64; nspec::Integer, reactions,
                         rate_of = nothing, rate_table = nothing,
                         nratec::Integer = 400,
                         temstart = 1.0, temend = 1.0e9) where {T}
    nreac = length(reactions)
    grid  = TempGrid(T; nratec = nratec, temstart = temstart, temend = temend)
    reac = zeros(Int, 3, nreac)
    prod = zeros(Int, 4, nreac)
    for (r, (rs, ps)) in enumerate(reactions)
        @assert length(rs) ≤ 3 "reaction $r has >3 reactants"
        @assert length(ps) ≤ 4 "reaction $r has >4 products"
        for (a, s) in enumerate(rs); reac[a, r] = s; end
        for (a, s) in enumerate(ps); prod[a, r] = s; end
    end
    # Rates: either a precomputed nratec×nreac table (the Krome path, tabulated
    # with invokelatest to dodge world-age), or built here from a `rate_of(r,T)`.
    if rate_table !== nothing
        @assert size(rate_table) == (nratec, nreac) "rate_table must be nratec×nreac"
        rate = convert(Matrix{T}, rate_table)
    else
        rate_of === nothing && error("provide either `rate_of` or `rate_table`")
        rate = Matrix{T}(undef, nratec, nreac)
        for r in 1:nreac, i in 1:nratec
            logttt = grid.logtem0 + (i - 1) * grid.dlogtem
            rate[i, r] = T(rate_of(r, exp(logttt)))
        end
    end
    return GenericNetwork{typeof(reac),typeof(rate),T}(nspec, nreac, reac, prod, rate, grid)
end

# column-indexed table interp (rate[:,r]) without allocating a view (GPU-safe)
@inline interp_col(rate, i::Int, r::Int, tdef::T) where {T} =
    @inbounds rate[i, r] + tdef * (rate[i + 1, r] - rate[i, r])

# per-cell rate of reaction r at column `c`, with the reactant product folded in
@inline function _reac_rate(n, c, reac, r, kr::T) where {T}
    @inbounds begin
        rate = kr
        s1 = reac[1, r]; s1 != 0 && (rate *= n[s1, c])
        s2 = reac[2, r]; s2 != 0 && (rate *= n[s2, c])
        s3 = reac[3, r]; s3 != 0 && (rate *= n[s3, c])
    end
    return rate
end

@kernel function _generic_step_kernel!(n, C, D, @Const(Tgas), reac, prod, rate,
                                       grid, dt, itmax::Int, nspec::Int, nreac::Int)
    c = @index(Global, Linear)          # one work-item per cell (column of n)
    T = eltype(n)
    @inbounds begin
        ii, tdef = table_index(log(max(Tgas[c], T(1))), grid)
        ttot = zero(T)
        for _ in 1:itmax
            for s in 1:nspec
                C[s, c] = zero(T); D[s, c] = zero(T)
            end
            # accumulate creation/destruction over reactions
            for r in 1:nreac
                kr = interp_col(rate, ii, r, tdef)
                R  = _reac_rate(n, c, reac, r, kr)        # full reaction rate
                # creation: add R to each product
                for a in 1:4
                    p = prod[a, r]; p != 0 && (C[p, c] += R)
                end
                # destruction coefficient: R / n_s for each reactant slot s
                for a in 1:3
                    s = reac[a, r]
                    if s != 0
                        ns = n[s, c]
                        D[s, c] += R / max(ns, T(RATE_TINY))
                    end
                end
            end
            # adaptive sub-step: 10% of the fastest fractional species change
            dtsub = dt - ttot
            for s in 1:nspec
                d = D[s, c]
                if d > T(RATE_TINY)
                    dtsub = min(dtsub, T(0.1) / d)
                end
            end
            dtsub = max(dtsub, T(1e-6) * dt)
            # semi-implicit BDF update of every species
            for s in 1:nspec
                n[s, c] = max((C[s, c] * dtsub + n[s, c]) / (one(T) + D[s, c] * dtsub),
                              T(RATE_TINY))
            end
            ttot += dtsub
            ttot ≥ dt * (one(T) - T(1e-6)) && break
        end
    end
end

"""
    allocate_workspace(be, T, nspec, ncells) -> (C, D)

Device scratch (`nspec×ncells` creation/destruction accumulators) for
[`generic_step!`](@ref). Allocate once and reuse across steps.
"""
function allocate_workspace(be, ::Type{T}, nspec::Integer, ncells::Integer) where {T}
    return device_zeros(be, T, (Int(nspec), Int(ncells))),
           device_zeros(be, T, (Int(nspec), Int(ncells)))
end

"""
    generic_step!(n, Tgas, net::GenericNetwork, C, D; dt, itmax=10000) -> n

Advance the `nspec×ncells` species matrix `n` by `dt` at the per-cell temperatures
`Tgas` (length `ncells`), using the mass-action network `net` and the scratch
`(C, D)` from [`allocate_workspace`](@ref). One work-item integrates one cell's
sub-cycled stiff network. Updates `n` in place.
"""
function generic_step!(n, Tgas, net::GenericNetwork, C, D; dt::Real, itmax::Integer = 10000)
    be = KA.get_backend(n)
    ncells = size(n, 2)
    _generic_step_kernel!(be)(n, C, D, Tgas, net.reac, net.prod, net.rate, net.grid,
                              dt, Int(itmax), net.nspec, net.nreac; ndrange = ncells)
    return n
end
