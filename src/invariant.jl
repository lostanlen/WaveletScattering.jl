abstract AbstractInvariant{T<:Number,D<:AbstractDomain}

immutable NullInvariant <: AbstractBank
end

immutable SumInvariant{T<:Number,D<:AbstractDomain} <: AbstractInvariant{T,D}
    domain::D
    pathkey::PathKey
    signaltype::Type{T}
end
SumInvariant(Ux::WaveletScattering.ScatteredBlob) =
    SumInvariant(Ux.nodes[collect(keys(Ux.nodes))[1]])

function SumInvariant{T<:Number}(node::AbstractNode{T})
    pathkey = node.ranges[1].first
    K = sum([r.first == pathkey for r in node.ranges])
    domain = SpatialDomain(K)
    signaltype = T
    SumInvariant(domain, pathkey, signaltype
end

function SumInvariant{T<:Number}(node::AbstractFourierNode{T})
    pathkey = node.ranges[1].first
    K = sum([r.first == pathkey for r in node.ranges])
    domain = FourierDomain(K)
    signaltype = T
    SumInvariant(domain, pathkey, signaltype)
end

function transform!{T<:Number}(
        destination::SubArray,
        invariant::SumInvariant{T,SpatialDomain},
        node::AbstractNode{T},
        dim::Int)
    destination[:] = sum(node.data, dim)
end
