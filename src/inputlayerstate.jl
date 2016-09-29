immutable InputLayerState <: AbstractScatteredLayerState
    blobs::Vector{Mocha.Blob}
    layer::InputLayer
end

function InputLayerState(backend::Mocha.Backend, layer::InputLayer)
    ranges = ntuple(k -> kthrange(layer, k), ndims(layer.data))
    blob = ScatteredBlob(
        DataStructures.SortedDict((Path() => Node(layer.data, ranges),)))
    blobs = Mocha.Blob[blob]
    return InputLayerState(blobs, layer)
end

function Mocha.setup(
        backend::Mocha.Backend,
        layer::InputLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    @assert length(diffs) == 0
    @assert length(inputs) == 0
    return InputLayerState(backend, layer)
end

kthrange(layer::InputLayer, k::Int) =
    (PathKey(layer.symbols[k]) => (1:1:size(layer.data, k)))
