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
        pairs = collect(inputs[idblob].nodes)
        blobs[idblob] = ScatteredBlob(Dict(map(layer.ρ, pairs)))
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
