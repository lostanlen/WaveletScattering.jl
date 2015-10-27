# WaveletLayer
# We adopt the same whitespace convention as in the Mocha code base
Mocha.@defstruct WaveletLayer Mocha.Layer (
    name :: AbstractString = "wavelets",
    (bottoms :: Vector{Symbol} = Symbol[], length(bottoms) > 0),
    (tops :: Vector{Symbol} = Symbol[], length(tops) == length(bottoms)),
)

Mocha.@characterize_layer(WaveletLayer,
    has_neuron => false,
    has_param => false,
    can_do_bp => true
)
