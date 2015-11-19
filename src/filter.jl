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

# product with scalar number is commutative
Base.(:*)(b::Number, ψ::AbstractFilter) = ψ * b

# element-wise multiplication operator ".*" with scalar falls back to "*"
Base.(:.*)(ψ::AbstractFilter, b::Number) = ψ * b

function spin!{T,D}(ψs::Array{AbstractFilter{T,D},3})
    (nΘs, nΧs, nJs) = size(ψs)
    idγs = range(1, nΘs, nΘ * nΧs * nJs)
    # Caution: partial linear indexing may be deprecated in the future
    pmap(idγ -> spin!(ψs[2:end, idγ], ψs[1, idγ]), idγs)
end
