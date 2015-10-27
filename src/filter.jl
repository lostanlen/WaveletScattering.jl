abstract AbstractFilter{T<:Number,K}

immutable Symmetric1DFilter{T<:Number} <: AbstractFilter{T,1}
    leg::Vector{T}
    zero::T
end
