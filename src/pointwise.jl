abstract AbstractPointwise

immutable Identity <: AbstractPointwise end

(ρ::Identity)(x::Any) = x

immutable Log1P{T<:AbstractFloat} <: AbstractPointwise
    threshold::T
    function (::Type{Log1P{T}})(threshold::T)
        (threshold<0) && throw(DomainError)
        new{T}(threshold)
    end
end

(ρ::Log1P)(blob_in::ScatteredBlob) = map(ρ, ScatteredBlob(blob_in.nodes))

function Base.abs{T,N}(nodes::DataStructures.SortedDict{
        Path,Node{T,N},Base.Order.ForwardOrdering})
    absnodes =
        DataStructures.SortedDict{Path,Node{T,N},Base.Order.ForwardOrdering}()
    for (path, nodevalue) in nodes
        absnodes[path] = Node{T,N}(abs(nodevalue.data), nodevalue.ranges)
    end
    return absnodes::DataStructures.SortedDict{
        Path,Node{T,N},Base.Order.ForwardOrdering}
end

Base.abs(Wx::ScatteredBlob) = ScatteredBlob(abs(Wx.nodes))

function Base.map(ρ::AbstractPointwise, innodes::SortedDict)
    outnodes =
        DataStructures.SortedDict{Path,Node,Base.Order.ForwardOrdering}()
    for path in keys(innodes)
        outnodes[path] = ρ(innodes[path])
    end
    return outnodes
end

function Base.map!(
        ρ::AbstractPointwise,
        blob_out::ScatteredBlob,
        blob_in::ScatteredBlob)
    @inbounds for id in eachindex(blob_in.nodes)
        map!(ρ, blob_out.nodes[id].data, blob_in.nodes[id].data)
    end
end
