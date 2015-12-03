# FourierLayer
# We adopt the same whitespace convention as in the Mocha code base
Mocha.@defstruct FourierLayer Mocha.Layer (
    name :: AbstractString = "fourier",
    (bottoms :: Vector{Symbol} = Symbol[], length(bottoms) > 0),
    (data :: Vector{Array} = Array[], length(data) == length(tops)),
)

Mocha.@characterize_layer(FourierLayer,
    has_neuron => false,
    has_param => false,
    can_do_bp => true
)
