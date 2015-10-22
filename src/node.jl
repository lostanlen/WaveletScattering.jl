# Node
abstract AbstractNode{T, N}

immutable FourierNode{T<:Complex,N} <: AbstractNode{T,N}
Base.complex{T<:Real}(::Type{T}) = Complex{T}
Base.complex{T<:Complex}(::Type{T}) = T

immutable FourierNode{T<:Number,C<:Complex, N} <: AbstractNode{T,N}
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

function FourierNode{T<:Number,N}(data::Array{T,N}, fourierdims,
                                  subscripts::NTuple{N, PathKey})
    ranges =
        ntuple(k -> PathRange(subscripts[k] => (1:size(data,k))), ndims(data))
    FourierNode(complex(data), fourierdims, ranges)
end

Base.fft!(node::FourierNode) = fft!(node.data, node.fourierdims)
Base.ifft!(node::FourierNode) = ifft!(node.data, node.fourierdims)
