"""A `Bank2D` is a two-dimensional wavelet filter bank, parametrized by
* `T`: numeric type of input, e.g. Float32, Float64.
* `D`: transform domain. Either `FourierDomain{2}` or `SpatialDomain{2}`.
* `G`: point group. Either `TrivialGroup` or `RotationGroup`.
* `W`: wavelet class, e.g. `Morlet` or `MexicanHat`.
Its fields are
* `ϕ`: low-pass filter, also called scaling function.
* `ψs`: 3d array of wavelets, indexed by spin, chroma, and octave.
* `behavior`: mutable behavior, e.g. oversampling.
* `spec`: immutable specifications, e.g. number of filters per octave.
To create a `Bank2D`
1. define a `Spec2D`,
2. if needed, provide behavior characteristics as keyword arguments.
Example:
spec = Spec2D(n_orientations = 8)
bank = Bank2D(spec)"""
immutable Bank2D{
        T<:Number,
        D<:PlaneDomains,
        G<:PlaneGroups,
        W<:RedundantWaveletClass} <: AbstractBank{T,D,G,W}
    ϕ::AbstractFilter{T,D}
    ψs::Array{AbstractFilter{T,D},3}
    behavior::Behavior{T}
    spec::Spec2D{T,D,G,W}

    function (::Type{Bank2D}){T,D,G,W}(
            spec::Spec2D{T,D,G,W},
            pathkey::PathKey = PathKey(:space) ;
            is_ϕ_applied::Bool = false,
            j_range::UnitRange{Int} = 0:(spec.n_octaves-1),
            log2_oversampling::Int = 0,
            max_log2_stride::Int = spec.n_octaves-1,
            weighting::AbstractWeighting = EqualWeighting())
        (nΘs, nΧs, nJs) = size(spec.ψmetas)
        ψs = Array(AbstractFilter{T,D}, (nΘs, nΧs, nJs))
        # TODO: write AbstractFilter{T}(y::Vector{T}, spec::AbstractSpec{T,FourierDomain{1}})
        #ψs[1, :, :] =
        #    pmap(AbstractFilter, spec.ψmetas[1, :, :], fill(spec, nΧs, nJs))
        #(nΘs > 1) && spin!(ψs)
        #ϕ = AbstractFilter(spec.ϕmeta, spec)
        #renormalize!(ϕ, ψs, spec)
        #behavior = Behavior(ϕ, ψs, spec, is_ϕ_applied, j_range,
        #    log2_oversampling, max_log2_stride, pathkey, weighting)
        #new{T,D,G,W}(ϕ, ψs, behavior, spec)
    end
end

Base.ndims(::Bank2D) = 2
