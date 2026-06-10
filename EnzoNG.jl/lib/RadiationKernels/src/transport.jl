# ── M1 moment transport: dimensionally-split global Lax-Friedrichs ────────────
# The (N, Fx, Fy, Fz) moment system is hyperbolic with max wavespeed c (the
# reduced speed of light). The global-Lax-Friedrichs (GLF) interface flux
#
#     F̂_{i+½} = ½(F_L + F_R) − ½ c (U_R − U_L)
#
# is the robust, positivity-friendly scheme used by RAMSES-RT. We sweep one axis
# at a time (Godunov split), each sweep reading the input field and writing a
# fresh output field (ping-pong) so neighbours are never half-updated -- the same
# discipline the hydro sweeps use. Boundaries are transmissive (neighbour index
# clamped to the cell itself).

"""
    RadField{A}

One photon group's M1 state on a grid: photon number density `N` and number-flux
components `Fx, Fy, Fz`, each a flat column-major `idim·jdim·kdim` array on one
backend. A multi-group field is a `Vector{RadField}`; the transport sweeps one
group at a time.
"""
struct RadField{A}
    N::A; Fx::A; Fy::A; Fz::A
end

RadField(be, ::Type{T}, dims::Dims) where {T} =
    RadField(device_zeros(be, T, dims), device_zeros(be, T, dims),
             device_zeros(be, T, dims), device_zeros(be, T, dims))

# GLF face flux between left (L) and right (R) states across an `axis`-face.
# Returns the flux 4-tuple (for N, Fx, Fy, Fz). `c` = (reduced) light speed.
@inline function _glf_face(NL::T, FxL::T, FyL::T, FzL::T,
                           NR::T, FxR::T, FyR::T, FzR::T, c::T, axis::Int) where {T}
    c2 = c * c
    # normal flux of N is the normal component of F
    FdL = axis == 1 ? FxL : axis == 2 ? FyL : FzL
    FdR = axis == 1 ? FxR : axis == 2 ? FyR : FzR
    fN  = T(0.5) * (FdL + FdR) - T(0.5) * c * (NR - NL)
    # flux of momentum component b across the axis-face is c²·P[axis,b]:
    # P[axis,axis] is the M1 diagonal, P[axis,b≠axis] the off-diagonal. Computed
    # component-by-component (no closures, so the GPU back-end inlines cleanly).
    PdiagL = pressure_tensor_diag(NL, FxL, FyL, FzL, c, axis)
    PdiagR = pressure_tensor_diag(NR, FxR, FyR, FzR, c, axis)
    PxL = axis == 1 ? PdiagL : pressure_tensor_offdiag(NL, FxL, FyL, FzL, c, axis, 1)
    PyL = axis == 2 ? PdiagL : pressure_tensor_offdiag(NL, FxL, FyL, FzL, c, axis, 2)
    PzL = axis == 3 ? PdiagL : pressure_tensor_offdiag(NL, FxL, FyL, FzL, c, axis, 3)
    PxR = axis == 1 ? PdiagR : pressure_tensor_offdiag(NR, FxR, FyR, FzR, c, axis, 1)
    PyR = axis == 2 ? PdiagR : pressure_tensor_offdiag(NR, FxR, FyR, FzR, c, axis, 2)
    PzR = axis == 3 ? PdiagR : pressure_tensor_offdiag(NR, FxR, FyR, FzR, c, axis, 3)
    fFx = T(0.5) * c2 * (PxL + PxR) - T(0.5) * c * (FxR - FxL)
    fFy = T(0.5) * c2 * (PyL + PyR) - T(0.5) * c * (FyR - FyL)
    fFz = T(0.5) * c2 * (PzL + PzR) - T(0.5) * c * (FzR - FzL)
    return fN, fFx, fFy, fFz
end

@kernel function _m1_sweep_kernel!(No, Fxo, Fyo, Fzo,
                                   @Const(N), @Const(Fx), @Const(Fy), @Const(Fz),
                                   c, dtdx, axis::Int, stride::Int,
                                   n_along::Int, idim::Int, jdim::Int, kdim::Int)
    gi, gj, gk = @index(Global, NTuple)
    T = eltype(N)
    idx = gi + (gj - 1) * idim + (gk - 1) * idim * jdim
    # position along the swept axis (1-based), to clamp the boundary faces
    pos = axis == 1 ? gi : axis == 2 ? gj : gk
    lm = pos == 1        ? idx : idx - stride       # left neighbour (transmissive)
    rp = pos == n_along  ? idx : idx + stride       # right neighbour (transmissive)
    @inbounds begin
        # left face: between lm and idx ; right face: between idx and rp
        fNL, fFxL, fFyL, fFzL = _glf_face(N[lm], Fx[lm], Fy[lm], Fz[lm],
                                          N[idx], Fx[idx], Fy[idx], Fz[idx], T(c), axis)
        fNR, fFxR, fFyR, fFzR = _glf_face(N[idx], Fx[idx], Fy[idx], Fz[idx],
                                          N[rp], Fx[rp], Fy[rp], Fz[rp], T(c), axis)
        No[idx]  = N[idx]  - T(dtdx) * (fNR  - fNL)
        Fxo[idx] = Fx[idx] - T(dtdx) * (fFxR - fFxL)
        Fyo[idx] = Fy[idx] - T(dtdx) * (fFyR - fFyL)
        Fzo[idx] = Fz[idx] - T(dtdx) * (fFzR - fFzL)
        # keep N non-negative (GLF is positivity-preserving under the CFL limit,
        # but f32 round-off near zero can dip slightly negative)
        if No[idx] < zero(T); No[idx] = zero(T); end
    end
end

"""
    m1_sweep!(out::RadField, in::RadField, c, dt, dx;
              axis, idim, jdim, kdim) -> out

One GLF transport sweep of `in` along `axis ∈ (1,2,3)` into `out`. `c` is the
(reduced) speed of light, `dx` the cell width; the CFL limit is `c·dt/dx ≤ 1/√D`
for a D-dimensional split (use `dt ≤ dx/(D·c)`). Boundaries are transmissive.
"""
function m1_sweep!(out::RadField, inp::RadField, c::Real, dt::Real, dx::Real;
                   axis::Integer, idim::Integer, jdim::Integer, kdim::Integer)
    be = KA.get_backend(inp.N)
    stride = axis == 1 ? 1 : axis == 2 ? Int(idim) : Int(idim) * Int(jdim)
    n_along = axis == 1 ? Int(idim) : axis == 2 ? Int(jdim) : Int(kdim)
    _m1_sweep_kernel!(be)(out.N, out.Fx, out.Fy, out.Fz, inp.N, inp.Fx, inp.Fy, inp.Fz,
                          c, dt / dx, Int(axis), stride, n_along,
                          Int(idim), Int(jdim), Int(kdim);
                          ndrange = (Int(idim), Int(jdim), Int(kdim)))
    return out
end

"""
    m1_step!(field::RadField, scratch::RadField, c, dt, dx;
             idim, jdim, kdim, ndim=3) -> field

Advance one group a full step by a Godunov split over the `ndim` axes, ping-ponging
between `field` and `scratch`. After an even number of sweeps the result is back in
`field` (returned); an odd `ndim` leaves it in `scratch`, which is copied back.
`dt` must satisfy the split CFL `dt ≤ dx/(ndim·c)`.
"""
function m1_step!(field::RadField, scratch::RadField, c::Real, dt::Real, dx::Real;
                  idim::Integer, jdim::Integer, kdim::Integer, ndim::Integer = 3)
    src, dst = field, scratch
    for axis in 1:ndim
        m1_sweep!(dst, src, c, dt, dx; axis = axis, idim = idim, jdim = jdim, kdim = kdim)
        src, dst = dst, src
    end
    if src !== field   # odd number of sweeps: result sits in scratch, copy home
        copyto!(field.N, src.N); copyto!(field.Fx, src.Fx)
        copyto!(field.Fy, src.Fy); copyto!(field.Fz, src.Fz)
    end
    return field
end
