# We adopt the same whitespace convention as in the Mocha code base
Mocha.@defstruct ScatteringLayer Mocha.Layer (
  name :: AbstractString = "scattering",
  (bottoms :: Vector{Symbol} = Symbol[], length(bottoms) > 0),
  (tops :: Vector{Symbol} = Symbol[], length(tops) == length(bottoms)),
  neuron :: ActivationFunction = Neurons.Identity(),
)
