abstract AbstractFilter{T<:Number,D<:AbstractDomain}

immutable Symmetric1DFilter{T<:Number} <: AbstractFilter{T,SpatialDomain{1}}
    leg::Vector{T}
    zero::T
end

immutable Asymmetric1DFilter{T<:Number} <: AbstractFilter{T,SpatialDomain{1}}
    neg::Vector{T}
    pos::Vector{T}
    zero::T
end

# product is commutative
Base.(:*)(b::Number, ψ::AbstractFilter) = ψ * b

# element-wise multiplication operator ".*" with scalar falls back to "*"
Base.(:.*)(ψ::AbstractFilter, b::Number) = ψ * b
