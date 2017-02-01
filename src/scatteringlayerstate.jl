# ScatteringLayerState
immutable ScatteringLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    blobs::Vector{B}
    layer::ScatteringLayer
end

function Mocha.setup(
        backend::Mocha.CPUBackend,
        layer::ScatteringLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    return ScatteringLayerState(backend, layer, inputs)
end
