abstract AbstractDomain
immutable FourierDomain{N} <: AbstractDomain
    dim::Val{N}
end
FourierDomain(dim::Int) = FourierDomain{dim}(Val{dim}())
immutable GraphDomain <: AbstractDomain end
immutable SpatialDomain{N} <: AbstractDomain
    dim::Val{N}
end
SpatialDomain(dim::Int) = SpatialDomain{dim}(Val{dim}())

typealias LineDomains Union{FourierDomain{1},SpatialDomain{1}}
typealias PlaneDomains Union{FourierDomain{2},SpatialDomain{2}}
