immutable FourierLayer <: Mocha.Layer
    bottoms::Vector{Symbol}
    flags::UInt32
    name::AbstractString
    pathkeys::Vector{PathKey}
    timelimit::Float64
    tops::Vector{Symbol}
end

function FourierLayer( ;
        bottoms::Vector{Symbol} = Symbol[],
        flags::UInt32 = FFTW.ESTIMATE,
        name::AbstractString = "fourier",
        pathkeys::Vector{PathKey} = PathKey[],
        timelimit::Float64 = Inf,
        tops::Vector{Symbol} = Symbol[])
    @assert length(bottoms) == 1
    @assert length(tops) == 1
    FourierLayer(bottoms, flags, name, pathkeys, timelimit, tops)
end

Mocha.can_do_bp(::FourierLayer) = false
Mocha.has_neuron(::FourierLayer) = false
Mocha.has_param(::FourierLayer) = false
