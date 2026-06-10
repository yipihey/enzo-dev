# ── Grackle Cloudy HDF5 table loader ──────────────────────────────────────────
# Reads a Grackle Cloudy cooling/heating table (CloudyData_*.h5 / *_metals.h5)
# into a ChemistryKernels.CloudyTable for the KA interpolation kernel. Grackle's
# convention: the LAST parameter axis is Temperature; Parameter1 is log10(n_H),
# Parameter2 (rank 3) the redshift; Cooling/Heating store log10 rates.
#
# HDF5 is row-major and HDF5.jl reverses dims into column-major Julia, so we
# re-orient the table to (n1[, nz], nT) and flatten row-major (n1-major, T
# fastest) to match interpolate_2D_g's `int_index = (i1-1)*n2 + i2`.

"""
    load_cloudy_table(path; with_heating=true, redshift=0.0) -> ChemistryKernels.CloudyTable

Load a Grackle Cloudy HDF5 table (`/Parameter1`=log10 n_H, `/Parameter2`=redshift
for rank 3, `/Temperature`=K, `/Cooling`,`/Heating`=log10 rates). Rank is inferred
from the `Cooling` dataset: 1 ⇒ log T only; 2 ⇒ (log n_H, log T); 3 ⇒ sliced at
the `redshift` plane nearest `redshift` to yield a rank-2 (log n_H, log T) table
(full 3-D z interpolation is a planned extension). `Temperature` → log10.
"""
function load_cloudy_table(path::AbstractString; with_heating::Bool = true,
                           redshift::Real = 0.0)
    isfile(path) || error("Cloudy table not found: $path")
    par1, par2, tgrid, cool, heat = h5open(path, "r") do f
        p1 = haskey(f, "Parameter1") ? Float64.(vec(read(f["Parameter1"]))) : Float64[]
        p2 = haskey(f, "Parameter2") ? Float64.(vec(read(f["Parameter2"]))) : Float64[]
        tg = Float64.(vec(read(f["Temperature"])))
        cl = Float64.(read(f["Cooling"]))
        ht = (with_heating && haskey(f, "Heating")) ? Float64.(read(f["Heating"])) : zeros(size(cl))
        (p1, p2, tg, cl, ht)
    end
    par2log = log10.(tgrid)            # last axis is always temperature
    nT = length(par2log)

    # Grackle stores the cooling/heating tables LINEARLY and log10's them on load
    # (cool1d_cloudy then exponentiates). Match that, flooring zeros (e.g. an
    # all-zero heating table) so 10^log → 0 rather than 10^0 = 1.
    log10floor(x) = log10(max(x, 1.0e-100))
    cool = log10floor.(cool); heat = log10floor.(heat)

    # rank-3 (n_H, z, T): slice at the redshift plane nearest `redshift`.
    # HDF5.jl reverses dims, so the Julia array is (nT, nz, n1) → slice dim 2.
    if ndims(cool) == 3
        iz = argmin(abs.(par2 .- redshift))
        cool = cool[:, iz, :]; heat = heat[:, iz, :]   # → (nT, n1)
    end

    if ndims(cool) == 1 || isempty(par1)
        # rank-1: temperature only
        n1 = nT
        dpar1 = (par2log[end] - par2log[1]) / (n1 - 1)
        return CloudyTable{Float64,Vector{Float64}}(1, n1, 1, 1,
            par2log, Float64[0.0], Float64[0.0], dpar1, 1.0, 1.0,
            Float64.(vec(cool)), Float64.(vec(heat)), with_heating)
    end

    # rank-2: orient to (n1 = density, nT = temperature)
    n1 = length(par1)
    C  = size(cool, 1) == n1 ? cool : permutedims(cool, (2, 1))   # → (n1, nT)
    H  = size(heat, 1) == n1 ? heat : permutedims(heat, (2, 1))
    @assert size(C) == (n1, nT) "unexpected Cooling shape $(size(cool)) vs ($n1,$nT)"
    flat(M) = vec(permutedims(M, (2, 1)))     # row-major: (i1-1)*nT + i2
    dpar1 = (par1[end] - par1[1]) / (n1 - 1)
    dpar2 = (par2log[end] - par2log[1]) / (nT - 1)
    return CloudyTable{Float64,Vector{Float64}}(2, n1, nT, 1,
        par1, par2log, Float64[0.0], dpar1, dpar2, 1.0,
        flat(C), flat(H), with_heating)
end
