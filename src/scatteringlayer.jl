immutable ScatteringLayer <: Mocha.Layer
    bank::AbstractBank
    bottoms::Vector{Symbol}
    name::AbstractString
    sibling_mask_factor::Float64
    tops::Vector{Symbol}
end

Mocha.can_do_bp(::ScatteringLayer) = false
Mocha.has_neuron(::ScatteringLayer) = false
Mocha.has_param(::ScatteringLayer) = false

function ScatteringLayer( ;
        bank::AbstractBank = NullBank(),
        bottoms::Vector{Symbol} = Symbol[],
        name::AbstractString = "scattering",
        sibling_mask_factor = 1.0,
        tops::Vector{Symbol} = Symbol[])
    @assert length(bottoms) == 2
    @assert length(tops) == 1
    ScatteringLayer(bank, bottoms, name, sibling_mask_factor, tops)
end
