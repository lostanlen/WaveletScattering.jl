# WaveletLayer
# We adopt the same whitespace convention as in the Mocha code base
Mocha.@defstruct WaveletLayer Mocha.Layer (
    name :: AbstractString = "wavelets",
    (bottoms :: Vector{Symbol} = Symbol[], length(bottoms) > 0),
    (tops :: Vector{Symbol} = Symbol[], length(tops) == length(bottoms)),
    neuron :: ActivationFunction = Mocha.Neurons.Identity()
)

Mocha.@characterize_layer(WaveletLayer,
    has_param => false,
    has_neuron => true,
    can_do_bp => true
)

# WaveletLayerState
immutable WaveletLayerState <: Mocha.LayerState
    bank::AbstractBank
    blobs::Vector{ScatteredBlob}
    blobs_diff::Vector{ScatteredBlob}
    layer::WaveletLayer
end

# ScatteredBlob
immutable ScatteredBlob{T<:Number, N} <: Mocha.Blob{T, N}
    nodes::Dict{Path, AbstractNode{T, N}}
    subscripts::NTuple{PathKey}
end

# Node
abstract AbstractNode{T, N}

end
