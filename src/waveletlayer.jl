immutable WaveletLayer <: Mocha.Layer
    bank::AbstractBank
    bottoms::Vector{Symbol}
    name::AbstractString
    tops::Vector{Symbol}
end

function WaveletLayer( ;
        bank::AbstractBank = NullBank(),
        bottoms::Vector{Symbol} = Symbol[],
        name::AbstractString = "wavelets",
        tops::Vector{Symbol} = Symbol[])
    @assert !isa(bank, NullBank)
    @assert length(bottoms) == 1
    @assert length(tops) == 1
    return WaveletLayer(bank, bottoms, name, tops)
end

Mocha.can_do_bp(::WaveletLayer) = false
Mocha.has_neuron(::WaveletLayer) = false
Mocha.has_param(::WaveletLayer) = false
