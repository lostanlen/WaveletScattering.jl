# ScatteredBlob
abstract AbstractScatteredBlob{T,N} <: Mocha.Blob{T,N}
abstract AbstractFourierBlob{T<:FFTW.fftwNumber,N} <: AbstractScatteredBlob{T,N}

immutable RealFourierBlob{T<:FFTW.fftwReal,N} <: AbstractFourierBlob{T,N}
    nodes::Dict{Path,RealFourierNode{T,N}}
    subscripts::NTuple{N,PathKey}
end

function RealFourierBlob{T<:FFTW.fftwReal,N}(node::RealFourierNode{T,N},
                                             subscripts::NTuple{N,PathKey})
    emptypath = Dict{PathKey,Int}()
    nodes = Dict(emptypath => node)
    RealFourierBlob{T,N}(nodes, subscripts)
end

Base.fft!(blob::AbstractFourierBlob) = pmap(fft!, values(blob.nodes))
Base.ifft!(blob::AbstractFourierBlob) = pmap(ifft!, values(blob.nodes))

function pathdepth(key::PathKey, refkey::PathKey)
    while ~isempty(key) && ~isempty(refkey) && (back(key) == back(refkey))
        pop!(key)
        pop!(refkey)
    end
    if isempty(refkey) && ~isempty(key)
        keyback = pop!(key)
        if isempty(key) && (keyback.symbol == :γ)
            return keyback.level
        end
    end
    return 0
end
function forward!(backend::Mocha.CPUBackend, blob_out::AbstractScatteredBlob,
                  bank::AbstractNonOrientedBank, blob_in::AbstractScatteredBlob)
    γkey = cons(Literal(:γ, 1), bank.behavior.pathkey)
    for γ in 0:(nGammas-1)
        ψ = bank.ψs[1 + γ]
        for (path_in, node_in) in input.nodes
            path_out = copy(path_in)
            path_out[γkey] = γ
            transform!(blob[path_out], blob_in[path_in], ψ)
        end
    end
end
