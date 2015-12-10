immutable InputLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    layer::InputLayer
    blob::B
end

function InputLayerState(backend::Backend, layer::InputLayer)
    ranges = ntuple(k -> kthrange(layer, k), ndims(layer.data))
    blob = ScatteredBlob(Path() => Node(layer.data, ranges))
    return InputLayerState(layer, blob)
end

kthrange(layer::InputLayer, k::Int) =
    (PathKey(layer.symbols[k]) => (1:1:size(layer.data, k)))
