# InputLayer
# We adopt the same whitespace convention as in the Mocha code base
Mocha.@defstruct InputLayer Mocha.Layer (
    name :: AbstractString = "scattering-input",
    top :: Symbol,
    data :: AbstractArray,
    (symbols :: Vector{Symbol}, length(symbols) == ndims(data))
)

Mocha.@characterize_layer(InputLayer,
    has_neuron => false,
    has_param => false,
    can_do_bp => true,
    is_source => true,
)
