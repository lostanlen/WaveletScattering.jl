# ScatteredBlob
immutable ScatteredBlob{T<:Number,N} <: Mocha.Blob{T, N}
    nodes::Dict{Path, AbstractNode{T,N}}
abstract AbstractScatteredBlob{T,N} <: Mocha.Blob{T,N}
    subscripts::NTuple{N,PathKey}
end

function ScatteredBlob{T<:Number,N}(x::AbstractArray{T,N},
                                    subscripts::NTuple{N,PathKey})
    emptypath = Dict{PathKey,Int}()
    nodes = Dict{Path, AbstractNode{T,N}}(emptypath => x)
    ScatteredBlob{T,N}(nodes, subscripts)
end

Base.fft!(blob::ScatteredBlob) = pmap(fft!, blob)

# WaveletLayer
# We adopt the same whitespace convention as in the Mocha code base
Mocha.@defstruct WaveletLayer Mocha.Layer (
    name :: AbstractString = "wavelets",
    (bottoms :: Vector{Symbol} = Symbol[], length(bottoms) > 0),
    (tops :: Vector{Symbol} = Symbol[], length(tops) == length(bottoms)),
    neuron :: ActivationFunction = Mocha.Neurons.Identity()
)

Mocha.@characterize_layer(WaveletLayer,
    has_param => false,
    has_neuron => true,
    can_do_bp => true
)

# WaveletLayerState
immutable WaveletLayerState <: Mocha.LayerState
    bank::AbstractBank
    blobs::Vector{ScatteredBlob}
    blobs_diff::Vector{ScatteredBlob}
    layer::WaveletLayer
end

function forward{B<:Mocha.Blob}(backend::Mocha.CPUBackend,
                                state::WaveletLayerState, inputs::Vector{B})
    for idblob in eachindex(inputs)
        forward!(state.blobs[idblob], backend, state, input[idblobd])
    end
end

function forward!(blob::ScatteredBlob, backend::CPUBackend,
                  state::WaveletLayerState, input::ScatteredBlob)

end
