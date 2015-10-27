# ScatteredBlob
abstract AbstractScatteredBlob{T,N} <: Mocha.Blob{T,N}
abstract AbstractFourierBlob{T,N} <: AbstractScatteredBlob{T,N}

immutable RealFourierBlob{T<:Real,N} <: AbstractFourierBlob{T,N}
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
