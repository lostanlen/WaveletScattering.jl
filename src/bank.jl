abstract AbstractBank{
        T<:Number,
        D<:AbstractDomain,
        G<:AbstractPointGroup,
        W<:RedundantWaveletClass}

immutable Bank1D{
        T<:Number,
        D<:LineDomains,
        G<:LineGroups,
        W<:RedundantWaveletClass} <: AbstractBank{T,D,G,W}
    ϕ::AbstractFilter{T,D}
    ψs::Array{AbstractFilter{T,D},2}
    behavior::Behavior
    spec::Spec1D{T,D,G,W}
    function call{T,D,G,W}(
            ::Type{Bank1D}, spec::Spec1D{T,D,G,W}, behavior::Behavior{G})
        ψs = pmap(AbstractFilter, metas, fill(spec, length(metas)))
        ψs = convert(Array{AbstractFilter{T,D,G,W},2}, ψs)
        ϕ = scalingfunction(spec)
        renormalize!(ϕ, ψs, metas, spec)
        behavior = Behavior(ϕ, ψs, spec,
            is_ϕ_applied, j_range, log2_oversampling, max_log2_stride)
        new{T,D,G,W}(ϕ, ψs, behavior, metas, spec)
    end
end
