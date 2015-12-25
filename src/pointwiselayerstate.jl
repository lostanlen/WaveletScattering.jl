immutable PointwiseLayerState{P<:AbstractPointwise}
    layer::PointwiseLayer{P}
    blobs::Vector{Mocha.Blob}
    blobs_diff::Vector{Mocha.Blob}
end

function PointwiseLayerState(
        backend::Mocha.CPUBackend,
        layer::PointwiseLayer,
        inputs::Vector{Mocha.Blob})
    blobs = Vector{Mocha.Blob}(length(inputs))
    for idblob in eachindex(inputs)
        pairs = collect(inputs[idblob].nodes)
        blobs[idblob] = ScatteredBlob(Dict(map(layer.ρ, pairs)))
    end
    return PointwiseLayerState(layer, blobs, ScatteredBlob[])
end

function forward(
        backend::Mocha.CPUBackend,
        layerstate::PointwiseLayerState,
        inputs::Vector{Mocha.Blob})
    for id in eachindex(inputs)
        map!(layerstate.ρ, layerstate.blobs[id], inputs[id])
    end
end
