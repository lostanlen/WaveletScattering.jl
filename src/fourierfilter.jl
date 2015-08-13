abstract AbstractFourierFilter{T<:Number}
abstract AbstractFourier1DFilter{T<:Number} <: AbstractFourierFilter

immutable Analytic1D{T<:Number}  <: AbstractFourier1DFilter{T}
    pos::Vector{T}
    posfirst::Int
end

immutable Coanalytic1D{T<:Number} <: AbstractFourier1DFilter{T}
    neg::Vector{T}
    neglast::Int
end

immutable Symmetric1D{T<:Number} <: AbstractFourier1DFilter{T}
    leg::Vector{T}
    zero::T
end

immutable Vanishing1D{T<:Number} <: AbstractFourier1DFilter{T}
    an::Analytic1D{T}
    coan::Coanalytic1D{T}
    midpoint::T
end
