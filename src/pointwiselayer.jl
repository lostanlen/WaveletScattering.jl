
immutable PointwiseLayer{P<:AbstractPointwise} <: Mocha.Layer
    name::AbstractString
    bottoms::Vector{Symbol}
    tops::Vector{Symbol}
    ρ::P
end

Mocha.can_do_bp(::PointwiseLayer) = true
Mocha.has_neuron(::PointwiseLayer) = false
Mocha.has_param(::PointwiseLayer) = false

function PointwiseLayer( ;
        name::AbstractString = "pointwise",
        bottoms::Vector{Symbol} = Symbol[],
        tops::Vector{Symbol} = Symbol[],
        ρ :: AbstractPointwise = Identity())
    @assert length(bottoms) > 0
    @assert length(tops) == length(bottoms)
    PointwiseLayer(name, bottoms, tops, ρ)
end
