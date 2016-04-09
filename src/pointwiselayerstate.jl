immutable PointwiseLayerState{P<:AbstractPointwise} <:
        AbstractScatteredLayerState
    layer::PointwiseLayer{P}
    blobs::Vector{Mocha.Blob}
    blobs_diff::Vector{Mocha.Blob}
end

function PointwiseLayerState(
        backend::Mocha.CPUBackend,
        layer::PointwiseLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    blobs = Vector{Mocha.Blob}(length(inputs))
    for idblob in eachindex(inputs)
        pairs = collect(inputs[idblob].nodes)
        blobs[idblob] = ScatteredBlob(
            DataStructures.SortedDict(map(layer.ρ, pairs)))
    end
    # TODO: build diffs
    return PointwiseLayerState(layer, blobs, diffs)
end

function Mocha.setup(
        backend::Mocha.Backend,
        layer::PointwiseLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    return PointwiseLayerState(backend, layer, inputs, diffs)
end

function forward(
        backend::Mocha.CPUBackend,
        layerstate::PointwiseLayerState,
        inputs::Vector{Mocha.Blob})
    for id in eachindex(inputs)
        map!(layerstate.ρ, layerstate.blobs[id], inputs[id])
    end
end
