# ── Cloudy metal-line cooling (KA port of Grackle's cool1d_cloudy_g.F) ────────
# Grackle's main addition over the Enzo primordial network is metal-line cooling
# interpolated from precomputed Cloudy tables. This is the GPU hot path: per cell
# a bi-/tri-linear interpolation in log space over (log n_H, [redshift], log T)
# of the Cooling and Heating tables, giving the metal contribution to edot.
#
# The interpolation is a verbatim port of Grackle's interpolate_1D_g/2D_g (log
# space, uniform parameter grids). The table is loaded by `Grackle.load_cloudy_table`
# (HDF5) and passed in as plain arrays so this library stays HDF5-free. Tables and
# kernels are precision-generic and run on CPU and Metal alike.

export CloudyTable, cloudy_metal_cooling!, cloudy_edot

const INV_LOG10 = 0.4342944819032518   # 1/ln(10)

"""
    CloudyTable{T,V}

A Grackle Cloudy cooling/heating table on a uniform log-parameter grid. `rank`
1 ⇒ T only; 2 ⇒ (log n_H, log T); 3 ⇒ (log n_H, redshift, log T). `par1/par2/par3`
are the grid axes (already in the interpolation's log units), `dpar*` their
uniform spacings, and `cooling`/`heating` the row-major flattened tables
(`gridDim1×…×gridDimRank`, last axis fastest — Grackle's layout). `with_heating`
toggles the heating term. Load with `Grackle.load_cloudy_table`.
"""
struct CloudyTable{T,V<:AbstractVector{T}}
    rank::Int
    n1::Int; n2::Int; n3::Int
    par1::V; par2::V; par3::V
    dpar1::T; dpar2::T; dpar3::T
    cooling::V; heating::V
    with_heating::Bool
end

# 1-D linear interpolation over a uniform-log grid (interpolate_1D_g)
@inline function _interp1d(x::T, par, dpar::T, n::Int, data, base::Int) where {T}
    @inbounds begin
        i = unsafe_trunc(Int, (x - par[1]) / dpar) + 1
        i = ifelse(i < 1, 1, ifelse(i > n - 1, n - 1, i))
        slope = (data[base + i + 1] - data[base + i]) / (par[i + 1] - par[i])
        return (x - par[i]) * slope + data[base + i]
    end
end

# 2-D bilinear interpolation (interpolate_2D_g): index1 over par1, index2 over par2,
# flattened data row-major as (index1-1)*n2 + index2.
@inline function _interp2d(x1::T, x2::T, table::CloudyTable{T}, data) where {T}
    @inbounds begin
        n2 = table.n2
        i1 = unsafe_trunc(Int, (x1 - table.par1[1]) / table.dpar1) + 1
        i1 = ifelse(i1 < 1, 1, ifelse(i1 > table.n1 - 1, table.n1 - 1, i1))
        i2 = unsafe_trunc(Int, (x2 - table.par2[1]) / table.dpar2) + 1
        i2 = ifelse(i2 < 1, 1, ifelse(i2 > n2 - 1, n2 - 1, i2))
        v = (zero(T), zero(T))
        # q=1,2 rows over par1
        vq1 = let idx = (i1 - 1) * n2 + i2
            s = (data[idx + 1] - data[idx]) / (table.par2[i2 + 1] - table.par2[i2])
            (x2 - table.par2[i2]) * s + data[idx]
        end
        vq2 = let idx = i1 * n2 + i2
            s = (data[idx + 1] - data[idx]) / (table.par2[i2 + 1] - table.par2[i2])
            (x2 - table.par2[i2]) * s + data[idx]
        end
        slope = (vq2 - vq1) / (table.par1[i1 + 1] - table.par1[i1])
        return (x1 - table.par1[i1]) * slope + vq1
    end
end

"""
    cloudy_edot(table, log_nH, log10T) -> edot_met

The metal cooling−heating rate (cooling-table units, *before* metallicity scaling)
at proper `log_nH = log10(n_H/cm³)` and `log10T`, by interpolating `table`.
`edot = -10^log_cool (+ 10^log_heat)` — the cool1d_cloudy assembly.
"""
@inline function cloudy_edot(table::CloudyTable{T}, log_nH::T, log10T::T) where {T}
    if table.rank == 1
        lc = _interp1d(log10T, table.par1, table.dpar1, table.n1, table.cooling, 0)
        e  = -exp10(lc)
        if table.with_heating
            lh = _interp1d(log10T, table.par1, table.dpar1, table.n1, table.heating, 0)
            e += exp10(lh)
        end
        return e
    else  # rank 2 (log n_H, log T)
        lc = _interp2d(log_nH, log10T, table, table.cooling)
        e  = -exp10(lc)
        if table.with_heating
            lh = _interp2d(log_nH, log10T, table, table.heating)
            e += exp10(lh)
        end
        return e
    end
end

# a 1-element placeholder table for the "no cloudy" path (passed but never read)
function _make_dummy_cloudy(::Type{T}) where {T}
    z = T[0]
    CloudyTable{T,Vector{T}}(1, 1, 1, 1, z, z, z, one(T), one(T), one(T), z, z, false)
end

@kernel function _cloudy_kernel!(edot_met, @Const(rhoH), @Const(logtem),
                                 @Const(metallicity), table, dom, zscale::Bool,
                                 i1::Int, j1::Int, idim::Int)
    gi, gj = @index(Global, NTuple)
    i = i1 + gi - 1; j = j1 + gj - 1
    idx = (j - 1) * idim + i
    T = eltype(edot_met)
    @inbounds begin
        log_nH = log10(max(rhoH[idx] * dom, T(RATE_TINY)))
        log10T = logtem[idx] * T(INV_LOG10)
        e = cloudy_edot(table, log_nH, log10T)
        edot_met[idx] = zscale ? e * metallicity[idx] : e
    end
end

"""
    cloudy_metal_cooling!(edot_met, rhoH, logtem, metallicity, table::CloudyTable;
                          dom, idim, i1, i2, j1=1, j2=1, zscale=true) -> edot_met

Fill `edot_met` with the per-cell metal cooling−heating rate (Grackle Cloudy
tables), the GPU hot path. `rhoH` is the H mass density (code units; ×`dom` →
proper n_H), `logtem` the natural-log gas temperature, `metallicity` the gas
metallicity (Z/Z_solar) when `zscale=true`. Mirrors `cool1d_cloudy_g.F`'s
`edot_met`; add it into the network's `edot`.
"""
function cloudy_metal_cooling!(edot_met, rhoH, logtem, metallicity, table::CloudyTable;
                               dom::Real, idim::Integer, i1::Integer, i2::Integer,
                               j1::Integer = 1, j2::Integer = 1, zscale::Bool = true)
    be = KA.get_backend(edot_met)
    ni = Int(i2 - i1 + 1); nj = Int(j2 - j1 + 1)
    _cloudy_kernel!(be)(edot_met, rhoH, logtem, metallicity, table, dom, zscale,
                        Int(i1), Int(j1), Int(idim); ndrange = (ni, nj))
    return edot_met
end
