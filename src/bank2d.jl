immutable Bank2D{
        T<:Number,
        D<:PlaneDomains,
        G<:PlaneGroups,
        W<:RedundantWaveletClass} <: AbstractBank{T,D,G,W}
    ϕ::AbstractFilter{T,D}
    ψs::Array{AbstractFilter{T,D},3}
    behavior::Behavior{T}
    spec::Spec2D{T,D,G,W}
end

Base.ndims(::Bank2D) = 2
