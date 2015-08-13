# Bank
abstract AbstractBank{T<:Number}
abstract AbstractNonOrientedBank{T<:Number} <: AbstractBank{T<:Number}
abstract AbstractOrientedBank{T<:Number} <: AbstractBank{T<:Number}

immutable FourierNonOriented1DBank{T<:Number} <: AbstractNonOrientedBank{T}
    ψs::Vector{AbstractFourier1DFilter{T}}
    ϕs::FourierSymmetric1DFilter{T}
    behavior::Behavior
    metas::Vector{NonOrientedMeta}
    spec::Abstract1DSpec{T}
end


# Behavior
type Behavior
    γ_range::UnitRange
    log2_oversampling::Int
    min_log2_resolution::Int
end
