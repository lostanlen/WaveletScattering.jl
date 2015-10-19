# ScatteredBlob
abstract AbstractScatteredBlob{T,N} <: Mocha.Blob{T,N}

immutable FourierScatteredBlob{T<:Number, N} <: AbstractScatteredBlob{T,N}
    nodes::Dict{Path,FourierNode{T,N}}
    subscripts::NTuple{N,PathKey}
end

function FourierScatteredBlob{T<:Number,N}(
        node::AbstractNode{T,N}, subscripts::NTuple{N,PathKey})
    emptypath = Dict{PathKey,Int}()
    nodes = Dict(emptypath => node)
    FourierScatteredBlob{T,N}(nodes, subscripts)
end

Base.fft!(blob::FourierScatteredBlob) = pmap(fft!, blob)

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
immutable WaveletLayerState{B<:AbstractScatteredBlob} <: Mocha.LayerState
    bank::AbstractBank
    blobs::Vector{B}
    blobs_diff::Vector{B}
    layer::WaveletLayer
end

function forward!(backend::Mocha.CPUBackend,
                  state::WaveletLayerState{FourierScatteredBlob},
                  inputs::Vector{FourierScatteredBlob})
    for idblob in eachindex(inputs)
        fft!(input[idblob])
        forward!(backend, state.blobs[idblob], state.bank, input[idblob])
        ifft!(state.blobs[idblob])
    end
end

function forward!(backend::Mocha.CPUBackend,
                  blob_out::AbstractScatteredBlob,
                  bank::AbstractNonOrientedBank,
                  blob_in::AbstractScatteredBlob)
    γkey = cons(Literal(:γ, 1), bank.behavior.pathkey)
    for γ in 0:(nGammas-1)
        ψ = bank.ψs[1+γ]
        for (path_in, node_in) in input.nodes
            path_out = copy(path_in)
            path_out[γkey] = γ
            transform!(blob[path_out], node_in[path_in], ψ)
        end
    end
end

end
