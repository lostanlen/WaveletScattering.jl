# Node
abstract AbstractNode{T, N}

immutable FourierNode{T<:Number,N} <: AbstractNode
    data::AbstractArray{T,N}
    data_ft::AbstractArray
    ranges::NTuple{PathRange, N}
end

function fft!(node::FourierNode, dims)
    blob.data_ft[:] = fft(blob.data, dims)
end

function fft!(blob::ScatteredBlob, dims)
    pmap(pair -> (pair.first, fft!(pair.second, dims)), blob)
end
