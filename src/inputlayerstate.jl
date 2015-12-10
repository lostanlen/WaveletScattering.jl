immutable InputLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    layer::InputLayer
    blob::B
end

function InputLayerState(backend::Backend, layer::InputLayer)
    ranges = ntuple(k -> kthrange(layer, k), ndims(layer.data))
    node = Node(layer.data, ranges)
end

kthrange(layer::InputLayer, k::Int) =
    (layer.symbols[k] => (1:1:size(layer.data, k)))s
