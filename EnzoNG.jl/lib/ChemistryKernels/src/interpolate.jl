# ── Log-temperature table lookup ──────────────────────────────────────────────
# The rate/cooling coefficients are tabulated on a uniform grid in log(T) over
# [temstart, temend] with `nratec` bins (calc_rates.F). cool1d_multi.F looks them
# up per cell with a clamped linear interpolation; these helpers are `@inline`
# device functions reused by every per-cell kernel, reproducing the exact
# `indixe`/`tdef` arithmetic Enzo uses (1-based, clamped to [1, nratec-1]).

"""
    TempGrid{T}

The uniform-in-log(T) tabulation grid shared by all coefficient tables:
`nratec` bins over `[temstart, temend]`, with `logtem0 = log(temstart)` and
`dlogtem = (log(temend)-log(temstart))/(nratec-1)`.
"""
struct TempGrid{T}
    nratec::Int
    logtem0::T
    dlogtem::T
    temstart::T
    temend::T
end

function TempGrid(::Type{T}; nratec::Integer, temstart, temend) where {T}
    l0 = log(T(temstart))
    dl = (log(T(temend)) - l0) / (nratec - 1)
    return TempGrid{T}(Int(nratec), l0, dl, T(temstart), T(temend))
end

"""
    table_index(logtem, grid) -> (i, tdef)

The 1-based lower bin index `i ∈ [1, nratec-1]` and linear interpolation weight
`tdef ∈ [0, 1]` for temperature `exp(logtem)`. Matches cool1d_multi.F:

    indixe = min(nratec-1, max(1, int((logtem - logtem0)/dlogtem) + 1))
    tdef   = (logtem - (logtem0 + (indixe-1)*dlogtem)) / dlogtem
"""
@inline function table_index(logtem::T, grid::TempGrid{T}) where {T}
    i = unsafe_trunc(Int, (logtem - grid.logtem0) / grid.dlogtem) + 1
    i = ifelse(i < 1, 1, ifelse(i > grid.nratec - 1, grid.nratec - 1, i))
    tlo = grid.logtem0 + T(i - 1) * grid.dlogtem
    tdef = (logtem - tlo) / grid.dlogtem
    return i, tdef
end

"`interp(tbl, i, tdef)` — linear interpolation `tbl[i] + tdef·(tbl[i+1]-tbl[i])`."
@inline interp(tbl, i::Int, tdef::T) where {T} =
    @inbounds tbl[i] + tdef * (tbl[i + 1] - tbl[i])
