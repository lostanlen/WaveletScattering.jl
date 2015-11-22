abstract AbstractPointwise
immutable Modulus <: AbstractPointwise end

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

function forward{B<:ScatteredBlob}(
        backend::Mocha.CPUBackend,
        layerstate::PointwiseLayerState,
        inputs::Vector{B})
    for id in eachindex(inputs)
        map!(layerstate.ρ, layerstate.blobs[id], inputs[id])
    end
end

map!(ρ::Modulus, blob_out::AbstractNode, blob_in::AbstractNode) =
    map!(abs, blob_out.data, blob_in.data)
