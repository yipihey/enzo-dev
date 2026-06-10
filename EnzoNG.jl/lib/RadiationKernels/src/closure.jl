# ── M1 (Levermore) closure ────────────────────────────────────────────────────
# Closes the moment system by expressing the radiation pressure tensor P in terms
# of the local (N, F). The reduced flux f = |F| / (c·N) ∈ [0, 1] interpolates
# between the isotropic (f→0, P = N/3·I) and free-streaming (f→1, P = N·n⊗n)
# limits. The Eddington factor is Levermore (1984):
#
#     χ(f) = (3 + 4 f²) / (5 + 2 √(4 − 3 f²))
#
# and  P = N · [ (1−χ)/2 · I + (3χ−1)/2 · n⊗n ],  n = F/|F|.

"""
    eddington_factor(f) -> χ

Levermore M1 Eddington factor for reduced flux `f = |F|/(cN) ∈ [0,1]`. Returns
`1/3` at `f=0` (isotropic) and `1` at `f=1` (free streaming). `f` is clamped to
`[0,1]` for safety against round-off pushing it slightly past the cone.
"""
@inline function eddington_factor(f::T) where {T}
    fc = ifelse(f < zero(T), zero(T), ifelse(f > one(T), one(T), f))
    return (T(3) + T(4) * fc * fc) / (T(5) + T(2) * sqrt(T(4) - T(3) * fc * fc))
end

"""
    pressure_tensor_diag(N, Fx, Fy, Fz, c, axis) -> P_axis,axis

The diagonal pressure-tensor component along `axis ∈ (1,2,3)` for the M1 closure,
i.e. `P[axis,axis]`. Used by the directional GLF flux for the normal-momentum
equation. `c` is the (reduced) speed of light setting `f = |F|/(cN)`.
"""
@inline function pressure_tensor_diag(N::T, Fx::T, Fy::T, Fz::T, c::T, axis::Int) where {T}
    Fmag2 = Fx*Fx + Fy*Fy + Fz*Fz
    Fmag  = sqrt(Fmag2)
    Nsafe = max(N, eps(T))
    f = Fmag / (c * Nsafe)
    χ = eddington_factor(f)
    na = axis == 1 ? Fx : axis == 2 ? Fy : Fz
    n2 = Fmag > eps(T) ? (na / Fmag)^2 : zero(T)
    return N * ((one(T) - χ) / T(2) + (T(3) * χ - one(T)) / T(2) * n2)
end

"""
    pressure_tensor_offdiag(N, Fx, Fy, Fz, c, a, b) -> P_a,b

Off-diagonal pressure component `P[a,b]` (`a ≠ b`) for the M1 closure, the flux of
the `b`-momentum across an `a`-face. `n_a n_b` direction cosines weight the
free-streaming part.
"""
@inline function pressure_tensor_offdiag(N::T, Fx::T, Fy::T, Fz::T, c::T, a::Int, b::Int) where {T}
    Fmag2 = Fx*Fx + Fy*Fy + Fz*Fz
    Fmag  = sqrt(Fmag2)
    Nsafe = max(N, eps(T))
    f = Fmag / (c * Nsafe)
    χ = eddington_factor(f)
    if Fmag ≤ eps(T)
        return zero(T)
    end
    na = a == 1 ? Fx : a == 2 ? Fy : Fz
    nb = b == 1 ? Fx : b == 2 ? Fy : Fz
    return N * (T(3) * χ - one(T)) / T(2) * (na * nb) / Fmag2
end
