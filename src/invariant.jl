abstract AbstractInvariant{T<:Number,D<:AbstractDomain}

immutable NullInvariant <: AbstractBank
end

immutable SumInvariant{T<:Number,D<:AbstractDomain} <: AbstractInvariant{T,D}
    domain::D
    pathkey::PathKey
    signaltype::Type{T}
end
SumInvariant(Ux::WaveletScattering.ScatteredBlob) =
    SumInvariant(Ux[collect(keys(Ux))[1]])

function SumInvariant{T<:Number}(node::AbstractNode{T})
    domain = SpatialDomain()
    pathkey = node.ranges[1].first
    signaltype = T
end

function SumInvariant{T<:Number}(node::AbstractFourierNode{T})
    domain = FourierDomain()
    pathkey = node.ranges[1].first
    signaltype = T
end

function transform!{T<:Number}(
        destination::SubArray,
        invariant::SumInvariant{T,SpatialDomain},
        node::AbstractNode{T,N},
        dim::Int)
    destination[:] = sum(node.data, dim)
end
