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
    SumInvariant(domain, pathkey, signaltype)
end

function SumInvariant{T<:Number}(node::AbstractFourierNode{T})
    pathkey = node.ranges[1].first
    K = sum([r.first == pathkey for r in node.ranges])
    domain = FourierDomain(K)
    signaltype = T
    SumInvariant(domain, pathkey, signaltype)
end

function (::Type{SumInvariant{T,SpatialDomain{K}}}){T,K}(
        nodes::DataStructures.SortedDict{
            WaveletScattering.Path,Node{T,N},Base.Order.ForwardOrdering})
    sumnodes =
        DataStructures.SortedDict{Path,Node{T,N},Base.Order.ForwardOrdering}()
    for (path, nodevalue) in nodes
        dims = find([r.first == invariant.pathkey for r in nodevalue.ranges])
        sumnodes[path] = Node{T,N}(sum(nodevalue.data, dims), nodevalue.ranges)
    end
    return sumnodes
end

call(invariant::SumInvariant, Ux::ScatteredBlob) = call(invariant, Ux.nodes)

function transform!{T<:Number}(
        destination::SubArray,
        invariant::SumInvariant{T,SpatialDomain},
        node::AbstractNode{T},
        dim::Int)
    destination[:] = sum(node.data, dim)
end
