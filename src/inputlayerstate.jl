immutable InputLayerState <: AbstractScatteredLayerState
    layer::InputLayer
    blobs::Vector{Mocha.Blob}
end

function InputLayerState(backend::Backend, layer::InputLayer)
    ranges = ntuple(k -> kthrange(layer, k), ndims(layer.data))
    blobs = [ ScatteredBlob(Dict(Path() => Node(layer.data, ranges))) ]
    return InputLayerState(layer, blobs)
end

function Mocha.setup(
        backend::Mocha.Backend,
        layer::InputLayer,
        inputs::Vector{Blob},
        diffs::Vector{Blob})
    return InputLayerState(backend, layer)
end

kthrange(layer::InputLayer, k::Int) =
    (PathKey(layer.symbols[k]) => (1:1:size(layer.data, k)))
