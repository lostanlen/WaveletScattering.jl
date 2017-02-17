abstract AbstractFilter{T<:Number,D<:AbstractDomain}

immutable Asymmetric1DFilter{T<:Number} <: AbstractFilter{T,SpatialDomain{1}}
    neg::Vector{T}
    pos::Vector{T}
    zero::T
end

immutable Symmetric1DFilter{T<:Number} <: AbstractFilter{T,SpatialDomain{1}}
    leg::Vector{T}
    zero::T
end

# product with scalar number is commutative
Base.:*(b::Number, ψ::AbstractFilter) = ψ * b

# element-wise multiplication operator ".*" with scalar falls back to "*"
Base.:.*(ψ::AbstractFilter, b::Number) = ψ * b

function spin!{T,D}(ψs::Array{AbstractFilter{T,D},3})
    (nΘs, nΧs, nJs) = size(ψs)
    for χ in 0:(nΧs-1), j in 0:(nJs-1)
        ψs[2:end, 1+χ, 1+j] = spin(ψs[1, 1+χ, 1+j])
    end
end
