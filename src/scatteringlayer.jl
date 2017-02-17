immutable ScatteringLayer <: Mocha.Layer
    bank::AbstractBank
    bottoms::Vector{Symbol}
    name::AbstractString
    previous_bank::AbstractBank
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
        previous_bank::AbstractBank = NullBank(),
        sibling_mask_factor = 1.0,
        tops::Vector{Symbol} = Symbol[])
    @assert length(bottoms) == 1
    @assert length(tops) == 1
    ScatteringLayer(bank, bottoms, name, previous_bank, sibling_mask_factor,
        tops)
end
