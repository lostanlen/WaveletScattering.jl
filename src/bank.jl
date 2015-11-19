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
    ψs::Array{AbstractFilter{T,D},3}
    behavior::Behavior
    spec::Spec1D{T,D,G,W}
    function call{T,D,G,W}(::Type{Bank1D},
            spec::Spec1D{T,D,G,W} ;
            is_ϕ_applied::Bool = false,
            j_range::UnitRange{Int} = 0:(spec.nOctaves-1),
            log2_oversampling::Int = 0,
            max_log2_stride::Int = spec.nOctaves-1)
        (nΘs, nΧs, nJs) = size(spec.ψmetas)
        ψs = Array(AbstractFilter{T,D}, (nΘs, nΧs, nJs))
        ψs[1, :, :] =
            pmap(AbstractFilter, spec.ψmetas[1, :, :], fill(spec, nΧs * nJs))
        (nΘs > 1) && spin!(ψs)
        ϕ = scalingfunction(spec)
        renormalize!(ϕ, ψs, spec)
        behavior = Behavior(ϕ, ψs, spec,
            is_ϕ_applied, j_range, log2_oversampling, max_log2_stride)
        new{T,D,G,W}(ϕ, ψs, behavior, spec)
    end
end
