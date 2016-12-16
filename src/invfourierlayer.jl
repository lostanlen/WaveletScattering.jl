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
    @assert length(bottoms) > 0
    @assert length(bottoms) == length(tops)
    InvFourierLayer(bottoms, flags, name, pathkeys, timelimit, tops)
end

Mocha.can_do_bp(::InvFourierLayer) = false
Mocha.has_neuron(::InvFourierLayer) = false
Mocha.has_param(::InvFourierLayer) = false
