immutable PointwiseLayer{P<:AbstractPointwise} <: Mocha.Layer
    ρ::P
    bottoms::Vector{Symbol}
    name::AbstractString
    tops::Vector{Symbol}
end

Mocha.can_do_bp(::PointwiseLayer) = false
Mocha.has_neuron(::PointwiseLayer) = false
Mocha.has_param(::PointwiseLayer) = false

function PointwiseLayer( ;
        ρ :: AbstractPointwise = Identity(),
        bottoms::Vector{Symbol} = Symbol[],
        name::AbstractString = "pointwise",
        tops::Vector{Symbol} = Symbol[])
    @assert length(bottoms) > 0
    @assert length(tops) == length(bottoms)
    PointwiseLayer(ρ, bottoms, name, tops)
end
