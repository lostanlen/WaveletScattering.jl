immutable FourierLayer <: Mocha.Layer
    name::AbstractString
    bottoms::Vector{Symbol}
    tops::Vector{Symbol}
    pathkeys::Vector{PathKey}
    flags::UInt32
    timelimit::Float64
end

function FourierLayer(
        name::AbstractString = "fourier",
        bottoms::Vector{Symbol} = Symbol[],
        tops::Vector{Symbol} = Symbol[],
        pathkeys::Vector{PathKey} = PathKey[],
        flags::UInt32 = FFTW.ESTIMATE,
        timelimit::Float64 = Inf)
    @assert length(bottoms) > 0
    @assert length(bottoms) == length(tops)
    FourierLayer(name, bottoms, tops, pathkeys, flags, timelimit)
end

Mocha.@characterize_layer(FourierLayer,
    has_neuron => false,
    has_param => false,
    can_do_bp => true
)
