# Node
abstract AbstractNode{T, N}

Base.complex{T<:Real}(::Type{T}) = Complex{T}
Base.complex{T<:Complex}(::Type{T}) = T

immutable FourierNode{T<:Number,C<:Complex,N} <: AbstractNode{T,N}
    data::Array{T,N}
    data_ft::Array{C,N}
    fourierdims::Vector{Int}
    ranges::NTuple{N, PathRange}
    function FourierNode(data::Array{T,N}, fourierdims, ranges)
        data_ft = complex(data)
        fft!(data_ft, fourierdims)
        new{T,complex(T),N}(data, data_ft, fourierdims, ranges)
    end
end

immutable Node{T<:Number,N} <: AbstractNode{T,N}
    data::Array{T,N}
    ranges::NTuple{N,PathRange}
end

function FourierNode{T<:Number,N}}(data::Array{T,N}, fourierdims, ranges)
    data_ft = complex(data)
    fft!(data_ft, fourierdims)
    FourierNode{T,complex(T),N}(data, data_ft, fourierdims, ranges)
end
function FourierNode{T<:Number,N}(data::Array{T,N}, fourierdims::Vector{Int},
                                  subscripts::NTuple{N, PathKey})
    ranges =
        ntuple(k -> PathRange(subscripts[k] => (1:1:size(data,k))), ndims(data))
    FourierNode(data, fourierdims, ranges)
end
FourierNode(data, fourierdims::Int, subscripts) =
    FourierNode(data, collect(fourierdims), subscripts)

function Base.fft!{T<:Real}(node::FourierNode{T})
    map!(complex, node.data_ft, node.data)
    fft!(node.data_ft, node.fourierdims)
end
function Base.fft!{T<:Complex}(node::FourierNode{T})
    copy!(node.data_ft, node.data)
    fft!(node.data_ft, node.fourierdims)
end
