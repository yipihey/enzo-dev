"""
    RadiationKernelsMetalExt

Package extension lighting up the Metal (Apple GPU) backend for `RadiationKernels`.
Loaded automatically when `Metal` is present; registers the `:metal` backend and
specialises the device-array helpers onto `MtlArray`. Metal is Float32-only.
"""
module RadiationKernelsMetalExt

using RadiationKernels
using Metal

function __init__()
    if Metal.functional()
        RadiationKernels.register_backend!(:metal, Metal.MetalBackend())
    end
end

RadiationKernels.device_zeros(::Metal.MetalBackend, ::Type{T}, dims::Dims) where {T} =
    Metal.zeros(T, dims)

end # module
