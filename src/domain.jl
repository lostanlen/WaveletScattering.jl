"""There are two subtypes of `AbstractDomain`:
* `FourierDomain{K}`: filtering is implemented as a product in the Fourier
domain. `K` is the dimension of the Fourier transform.
* `SpatialDomain`: filtering is implemented as a convolution in the spatial
domain, or as a wavelet lifting scheme when it is possible. `K` is the number of
subscripts over which the convolution operates.
"""
abstract AbstractDomain

immutable FourierDomain{K<:Int} <: AbstractDomain
    dim::Val{K}
    function FourierDomain{D<:Int}(dim::Val{D})
        @assert (D > 0)
        new(FourierDomain{D}(dim))
    end
end
FourierDomain(dim::Int) = FourierDomain{dim}(Val{dim}())


immutable SpatialDomain{K<:Int} <: AbstractDomain
    dim::Val{K}
    function SpatialDomain{D<:Int}(dim::Val{D})
        @assert (D > 0)
        new(SpatialDomain{D}(dim))
    end
end
SpatialDomain(dim::Int) = SpatialDomain{dim}(Val{dim}())

typealias LineDomains Union{FourierDomain{1},SpatialDomain{1}}
typealias PlaneDomains Union{FourierDomain{2},SpatialDomain{2}}
