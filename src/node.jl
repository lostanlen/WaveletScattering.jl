# Node
abstract AbstractNode{T, N}

immutable FourierNode{T<:Number,N} <: AbstractNode
    data::AbstractArray{T,N}
    data_ft::AbstractArray
    ranges::NTuple{N, PathRange}
end

function fft!(node::FourierNode, dims)
    node.data_ft[:] = fft(node.data, dims)
end
