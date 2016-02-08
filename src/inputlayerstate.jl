immutable InputLayerState <: AbstractScatteredLayerState
    blobs::Vector{Mocha.Blob}
    layer::InputLayer
end

function InputLayerState(backend::Backend, layer::InputLayer)
    ranges = ntuple(k -> kthrange(layer, k), ndims(layer.data))
    blob = ScatteredBlob(Dict(Path() => Node(layer.data, ranges)))
    blobs = Mocha.Blob[blob]
    return InputLayerState(blobs, layer)
end

function Mocha.setup(
        backend::Mocha.Backend,
        diffs::Vector{Blob},
        inputs::Vector{Blob},
        layer::InputLayer)
    @assert length(inputs) == 0
    return InputLayerState(backend, layer)
end

kthrange(layer::InputLayer, k::Int) =
    (PathKey(layer.symbols[k]) => (1:1:size(layer.data, k)))
kthrange(syms::Vector{Symbol}, x::AbstractArray, k::Int) =
    (PathKey(syms[k]) => (1:1:size(x, k)))
