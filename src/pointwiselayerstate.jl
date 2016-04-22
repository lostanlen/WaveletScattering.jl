immutable PointwiseLayerState{P<:AbstractPointwise} <:
        AbstractScatteredLayerState
    blobs::Vector{Mocha.Blob}
    layer::PointwiseLayer{P}
end

function PointwiseLayerState(
        backend::Mocha.CPUBackend,
        layer::PointwiseLayer,
        inputs::Vector{Mocha.Blob})
    blobs = Vector{Mocha.Blob}(length(inputs))
    for idblob in eachindex(inputs)
        pairs = collect(inputs[idblob].nodes)
        blobs[idblob] = ScatteredBlob(
            DataStructures.SortedDict(map(layer.ρ, pairs)))
    end
    return PointwiseLayerState(blobs, layer)
end

function Mocha.setup(
        backend::Mocha.CPUBackend,
        layer::PointwiseLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    @assert length(inputs) == length(diffs)
    return PointwiseLayerState(backend, layer, inputs)
end

function forward(
        backend::Mocha.CPUBackend,
        layerstate::PointwiseLayerState,
        inputs::Vector{Mocha.Blob})
    for id in eachindex(inputs)
        map!(layerstate.ρ, layerstate.blobs[id], inputs[id])
    end
end
