# Node
abstract AbstractNode{T, N}

immutable InplaceFourierNode{T<:Number,N} <: AbstractNode
    data::Array{T,N}
    fourierdims
    ranges::NTuple{N, PathRange}
end

function InplaceFourierNode{T<:Complex,N}(data::Array{T,N},
                                         fourierdims,
                                         subscripts::NTuple{N, PathKey})
    ranges = ntuple(k -> PathRange(subscripts(k),1:size(data,k)), ndims(data))
    FourierNode(data, fourierdims, ranges)
end

fft!(node::InplaceFourierNode) = fft(node.data, node.fourierdims)
