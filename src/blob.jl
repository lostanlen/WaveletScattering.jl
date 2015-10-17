# ScatteredBlob
immutable ScatteredBlob{T<:Number,N} <: Mocha.Blob{T, N}
    nodes::Dict{Path, AbstractNode{T,N}}
    subscripts::NTuple{N,PathKey}
end

function ScatteredBlob{T<:Number,N}(x::AbstractArray{T,N},
                                    subscripts::NTuple{N,PathKey})
    emptypath = Dict{PathKey,Int}()
    nodes = Dict{Path, AbstractNode{T,N}}(emptypath => x)
    ScatteredBlob{T,N}(nodes, subscripts)
end

fft!(blob::ScatteredBlob) = pmap(fft!, blob)

function forward!(blob::ScatteredBlob, backend::CPUBackend,
                  state::WaveletLayerState, input::ScatteredBlob)

end
