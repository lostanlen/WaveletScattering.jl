abstract AbstractFilter{T<:Number}

immutable Symmetric1DFilter{T<:Number} <: AbstractFilter{T}
    leg::Vector{T}
    zero::T
end
