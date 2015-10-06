# WaveletLayer
# We adopt the same whitespace convention as in the Mocha code base
Mocha.@defstruct WaveletLayer Mocha.Layer (
    name :: AbstractString = "scattering",
    (bottoms :: Vector{Symbol} = Symbol[], length(bottoms) > 0),
    (tops :: Vector{Symbol} = Symbol[], length(tops) == length(bottoms)),
    neuron :: ActivationFunction = Mocha.Neurons.Identity(),
)

Mocha.@characterize_layer(WaveletLayer,
    has_param => false,
    has_neuron => true,
    can_do_bp => true
)

# WaveletLayerState
immutable WaveletLayerState{B<:AbstractBank} <: Mocha.LayerState
    layer::WaveletLayer
    blobs::Vector{Mocha.Blob}
    blobs_diff::Vector{Mocha.Blob}
    bank::B
# AbstractScatteredBlob
abstract AbstractScatteredBlob{T} <: Mocha.Blob{T}
end
