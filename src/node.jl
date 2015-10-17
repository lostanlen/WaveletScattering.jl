# Node
abstract AbstractNode{T, N}

immutable InplaceFourierNode{T<:Complex,N} <: AbstractNode
    data::Array{T,N}
    fourierdims
    ranges::NTuple{N, PathRange}
end

function InplaceFourierNode{T<:Number,N}(data::Array{T,N},
                                         fourierdims,
                                         subscripts::NTuple{N, PathKey})
    ranges =
        ntuple(k -> PathRange(subscripts[k] => (1:size(data,k))), ndims(data))
    InplaceFourierNode(complex(data), fourierdims, ranges)
end

fft!(node::InplaceFourierNode) = fft!(node.data, node.fourierdims)
