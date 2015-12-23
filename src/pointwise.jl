abstract AbstractPointwise

function call{T,N}(ρ::AbstractPointwise, node::AbstractNode{T,N})
    return Node(ρ(node.data), node.ranges)
end

immutable Identity <: AbstractPointwise
end

call{T,N}(ρ::Identity, data::AbstractArray{T,N}) = data

immutable Modulus <: AbstractPointwise
end

call{T,N}(ρ::Modulus, data::AbstractArray{T,N}) = abs(data)

immutable SquaredModulus <: AbstractPointwise
end

call{T,N}(ρ::SquaredModulus, data::AbstractArray{T,N}) = abs2(data)

immutable Log1P{T<:AbstractFloat} <: AbstractPointwise
    threshold::T
end

call{T,N}(ρ::Log1P, data::AbstractArray{T,N}) = log1p(ρ.threshold * data)

immutable PointwiseLayer{P<:AbstractPointwise} <: Mocha.Layer
    name::AbstractString
    bottoms::Vector{Symbol}
    tops::Vector{Symbol}
    ρ::P
end

function PointwiseLayer( ;
        name::AbstractString = "pointwise",
        bottoms::Vector{Symbol} = Symbol[],
        tops::Vector{Symbol} = Symbol[],
        ρ :: AbstractPointwise = Identity())
    @assert length(bottoms) > 0
    @assert length(tops) == length(bottoms)
    PointwiseLayer(name, bottoms, tops, ρ)
end

immutable PointwiseLayerState{BLOB<:ScatteredBlob,P<:AbstractPointwise}
    layer::PointwiseLayer{P}
    blobs::Vector{BLOB}
    blobs_diff::Vector{BLOB}
end

function PointwiseLayerState{B<:ScatteredBlob}(
        backend::Mocha.CPUBackend,
        layer::PointwiseLayer,
        inputs::Vector{B})
    blobs = Vector{ScatteredBlob}(length(inputs))
    for idblob in eachindex(inputs)
        blobs[idblob] = ScatteredBlob(map(layer.ρ, inputs[idblob].nodes))
    end
    return PointwiseLayerState(layer, blobs, ScatteredBlob[])
end

function forward{B<:ScatteredBlob}(
        backend::Mocha.CPUBackend,
        layerstate::PointwiseLayerState,
        inputs::Vector{B})
    for id in eachindex(inputs)
        map!(layerstate.ρ, layerstate.blobs[id], inputs[id])
    end
end

function map!(
        ρ::AbstractPointwise,
        blob_out::ScatteredBlob,
        blob_in::ScatteredBlob)
    @inbounds for id in eachindex(blob_in.nodes)
        map!(ρ, blob_out.nodes[id].data, blob_in.nodes[id].data)
    end
end

map!{T<:Real}(ρ::Modulus, data_out::Array{T}, data_in::Array{Complex{T}}) =
    map!(abs, data_out, data_in)

function map!{T<:Real}(ρ::Log1P{T}, data_out::Array{T}, data_in::Array{T})
    @inbounds @fastmath for id in eachindex(data_in)
        data_out[id] = data_in[id] * ρ.threshold
        data_out[id] = log1p(data_out[id])
    end
end
