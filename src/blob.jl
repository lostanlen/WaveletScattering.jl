# ScatteredBlob
abstract AbstractScatteredBlob{T,D<:AbstractDomain,N} <: Mocha.Blob{T,N}

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

pathdepth(path::Path, refkey::PathKey) =
    mapreduce(path -> pathdepth(path, refkey), max, 1, keys(path))

function pathdepth(key::PathKey, refkey::PathKey)
    while ~isempty(key) && ~isempty(refkey) && (back(key) == back(refkey))
        pop!(key)
        pop!(refkey)
    end
    if isempty(refkey) && ~isempty(key)
        keyback = pop!(key)
        isempty(key) && (keyback.symbol == :γ) && return (1 + keyback.level)
    end
    return 1
end

function pathdepth(blob::AbstractFourierBlob, refkey::PathKey)
    anypath = keys(blob.nodes)[1]
    pathdepth_dictlevel = pathdepth(anypath)
    pathdepth_tensorlevel = pathdepth(blob.nodes[anypath])
    return max(pathdepth_dictlevel, pathdepth_tensorlevel)
end

function forward!(backend::Mocha.CPUBackend, blob_out::AbstractScatteredBlob,
                  bank::AbstractNonOrientedBank, blob_in::AbstractScatteredBlob)
    pathdepth(blob_in, bank.behavior.pathkey)

    map(node -> pathdepth(bank.behavior.pathkey, keys(blob_in.nodes))
    γkey = cons(Literal(:γ, 1), bank.behavior.pathkey)
    for j in bank.behavior.j_range
        for χ in 0:(bank.spec.nFilters_per_octave-1)
            γ = j * nFilters_per_octave + χ
            ψ = bank.ψs[1 + γ]
            for (path_in, node_in) in input.nodes
                path_out = copy(path_in)
                path_out[γkey] = γ
                transform!(blob[path_out], blob_in[path_in], ψ)
            end
        end
    end
end
