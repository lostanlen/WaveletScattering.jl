# Node
abstract AbstractNode{T, N}

Base.complex{T<:Real}(::Type{T}) = Complex{T}
Base.complex{T<:Complex}(::Type{T}) = T

immutable FourierNode{T<:Number,C<:Complex,N} <: AbstractNode{T,N}
    data::Array{T,N}
    data_ft::Array{C,N}
    fourierdims
    ranges::NTuple{N, PathRange}
    fourier FourierNode(data::Array{T,N}, fourierdims, ranges)
        data_ft = complex(data)
        fft!(data_ft, fourierdims)
        new{T, complex(T), N}(data, data_ft, fourierdims, ranges)
    end
end

immutable Node{T<:Number,N} <: AbstractNode{T,N}
    data::Array{T,N}
    ranges::NTUple{N,PathRange}
end

function FourierNode{T<:Number,N}(data::Array{T,N}, fourierdims,
                                  subscripts::NTuple{N, PathKey})
    ranges =
        ntuple(k -> PathRange(subscripts[k] => (1:size(data,k))), ndims(data))
    FourierNode(complex(data), fourierdims, ranges)
end

function Base.fft!(node::FourierNode{T<:Real})
    map!(complex, node.data_ft, node.data)
    fft!(node.data_ft, node.fourierdims)
end
function Base.fft!(node::FourierNode{T<:Complex})
    copy!(node.data_ft, node.data)
    fft!(node.data_ft, node.fourierdims)
end
