immutable InvFourierLayer <: Mocha.Layer
    bottoms::Vector{Symbol}
    flags::UInt32
    name::AbstractString
    pathkeys::Vector{PathKey}
    timelimit::Float64
    tops::Vector{Symbol}
end

function InvFourierLayer( ;
        bottoms::Vector{Symbol} = Symbol[],
        flags::UInt32 = FFTW.ESTIMATE,
        name::AbstractString = "invfourier",
        pathkeys::Vector{PathKey} = PathKey[],
        timelimit::Float64 = Inf,
        tops::Vector{Symbol} = Symbol[])
    @assert length(bottoms) == 1
    @assert length(tops) == 1
    InvFourierLayer(bottoms, flags, name, pathkeys, timelimit, tops)
end

Mocha.can_do_bp(::InvFourierLayer) = false
Mocha.has_neuron(::InvFourierLayer) = false
Mocha.has_param(::InvFourierLayer) = false
