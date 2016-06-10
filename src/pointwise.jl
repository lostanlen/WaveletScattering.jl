abstract AbstractPointwise

immutable Identity <: AbstractPointwise end

call{T,N}(ρ::Identity, data::AbstractArray{T,N}) = data

immutable Log1P{T<:AbstractFloat} <: AbstractPointwise
    threshold::T
    function call{T}(::Type{Log1P}, threshold::T)
        (threshold<0) && throw(DomainError)
        new{T}(threshold)
    end
end

function Base.abs{NODE<:AbstractNode}(nodes::DataStructures.SortedDict{
        WaveletScattering.Path,NODE,Base.Order.ForwardOrdering})
    nodevalues = collect(values(nodes))
    isempty(nodevalues) && return nodes
    T = real(eltype(nodevalues[1].data))
    dim = ndims(nodevalues[1].data)
    absnodes =
        DataStructures.SortedDict{Path,Node{T,dim},Base.Order.ForwardOrdering}()
    for (path, nodevalue) in nodes
        absnodes[path] = Node{T,dim}(abs(nodevalue.data), nodevalue.ranges)
    end
end

Base.abs(Wx::ScatteredBlob) = ScatteredBlob(abs(Wx.nodes))

function Base.call(ρ::AbstractPointwise, blob_in::ScatteredBlob)
    blob_out = ScatteredBlob(map(ρ, blob_in.nodes))
end

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
