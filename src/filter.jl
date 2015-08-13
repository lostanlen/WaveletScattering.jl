abstract AbstractFilter{T<:Number}

immutable Symmetric1DFilter{T<:Number} <: AbstractFilter
    leg::Vector{T}
    zero::T
end
