# InputLayer
# We adopt the same whitespace convention as in the Mocha code base
Mocha.@defstruct InputLayer Mocha.Layer (
    data :: AbstractArray = [],
    name :: AbstractString = "signal",
    (symbols :: Vector{Symbol} = Symbol[], length(symbols) == ndims(data)),
    tops :: Vector{Symbol} = :data,
)

Mocha.can_do_bp(::InputLayer) = false
Mocha.has_neuron(::InputLayer) = false
Mocha.has_param(::InputLayer) = false
Mocha.is_source(::InputLayer) = true
