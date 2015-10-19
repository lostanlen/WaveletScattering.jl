# Node
abstract AbstractNode{T, N}

immutable FourierNode{T<:Complex,N} <: AbstractNode
    data::Array{T,N}
    fourierdims
    ranges::NTuple{N, PathRange}
end

function FourierNode{T<:Number,N}(data::Array{T,N}, fourierdims,
                                  subscripts::NTuple{N, PathKey})
    ranges =
        ntuple(k -> PathRange(subscripts[k] => (1:size(data,k))), ndims(data))
    FourierNode(complex(data), fourierdims, ranges)
end

Base.fft!(node::FourierNode) = fft!(node.data, node.fourierdims)
