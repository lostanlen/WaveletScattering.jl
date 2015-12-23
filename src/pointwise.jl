abstract AbstractPointwise

function call{T,N}(ρ::AbstractPointwise, node::AbstractNode{T,N})
    return Node(ρ(node.data), node.ranges)
end

immutable Modulus <: AbstractPointwise
end

call{T,N}(ρ::Modulus, data::AbstractArray{T,N}) = abs(data)

immutable Log1P{T<:AbstractFloat} <: AbstractPointwise
    threshold::T
end

Mocha.@defstruct PointwiseLayer Mocha.Layer (
    name :: AbstractString = "pointwise",
    (bottoms :: Vector{Symbol} = Symbol[], length(bottoms) > 0),
    (tops :: Vector{Symbol} = Symbol[], length(tops) == length(bottoms)),
)

immutable PointwiseLayerState{BLOB<:ScatteredBlob,P<:AbstractPointwise}
    layer::PointwiseLayer
    blobs::Vector{BLOB}
    blobs_diff::Any
    ρ::P
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
