# InputLayer
# We adopt the same whitespace convention as in the Mocha code base
Mocha.@defstruct InputLayer Mocha.Layer (
    name :: AbstractString = "signal",
    tops :: Vector{Symbol} = :data,
    data :: AbstractArray = [],
    (symbols :: Vector{Symbol} = Symbol[], length(symbols) == ndims(data))
)

Mocha.can_do_bp(::InputLayer) = true
Mocha.has_neuron(::InputLayer) = false
Mocha.has_param(::InputLayer) = false
Mocha.is_source(::InputLayer) = true
