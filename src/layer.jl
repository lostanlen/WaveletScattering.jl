# ScatteredBlob
abstract AbstractScatteredBlob{T,N} <: Mocha.Blob{T,N}
abstract AbstractFourierBlob{T,N} <: AbstractScatteredBlob{T,N}

immutable RealFourierBlob{T<:Real,N} <: AbstractScatteredBlob{T,N}
    nodes::Dict{Path,RealFourierNode{T,N}}
    subscripts::NTuple{N,PathKey}
end

abstract AbstractPointwise
immutable Modulus <: AbstractPointwise end

map!(ρ::Modulus, blob_in::AbstractNode, blob_out::AbstractNode) =
    map!(abs, blob_in.data, blob_out.data)

function RealFourierBlob{T<:Number,N}(node::RealFourierNode{T,N},
                                      subscripts::NTuple{N,PathKey})
    emptypath = Dict{PathKey,Int}()
    nodes = Dict(emptypath => node)
    RealFourierBlob{T,N}(nodes, subscripts)
end

fft!(blob::AbstractFourierBlob) = pmap(fft!, values(blob.nodes))

# WaveletLayer
# We adopt the same whitespace convention as in the Mocha code base
Mocha.@defstruct WaveletLayer Mocha.Layer (
    name :: AbstractString = "wavelets",
    (bottoms :: Vector{Symbol} = Symbol[], length(bottoms) > 0),
    (tops :: Vector{Symbol} = Symbol[], length(tops) == length(bottoms)),
    neuron :: ActivationFunction = Mocha.Neurons.Identity()
)

Mocha.@characterize_layer(WaveletLayer,
    has_neuron => false,
    has_param => false,
    can_do_bp => true
)

# WaveletLayerState
immutable WaveletLayerState{B<:AbstractScatteredBlob} <: Mocha.LayerState
    bank::AbstractBank
    blobs::Vector{B}
    blobs_diff::Vector{B}
    layer::WaveletLayer
end

function forward!(backend::Mocha.CPUBackend, state::WaveletLayerState,
                  ρ::AbstractPointwise, inputs::Vector)
    @inbounds for idblob in eachindex(inputs)
        map!(ρ, state.blobs[idblob], inputs[idblob])
    end
end

function forward!(backend::Mocha.CPUBackend,
                  state::WaveletLayerState{AbstractFourierBlob},
                  inputs::Vector{AbstractFourierBlob})
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

function transform!(node_in::AbstractFourierNode,
                    node_out::AbstractFourierNode,
                    ψ::FullResolution1DFilter)
    inds = fill(Colon(), ndims(node_in))
    N = length(ψ.coeff)
    for ω in 0:(N-1)
        inds[node_in.fourierdims[1]] = 1 + ω
        view_in = view(node_in, inds)
        view_out = view(node_out, inds)
        for id in eachindex(view_in)
            view_out[id] = view_in[id] * ψ[1 + ω]
        end
    end
end
