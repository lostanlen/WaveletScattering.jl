abstract AbstractFilter{T<:Number,N} <: AbstractArray{T,N}

immutable Symmetric1DFilter{T<:Number} <: AbstractFilter{T,1}
    leg::Vector{T}
    zero::T
end
