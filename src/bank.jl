abstract AbstractBank{
    T<:Number,
    D<:AbstractDomain,
    G<:AbstractPointGroup,
    W<:RedundantWaveletClass}

"""A `Bank1D` is a one-dimensional wavelet filter bank, parametrized by
* `T`: numeric type of input, e.g. Float32, Float64.
* `D`: transform domain. Either `FourierDomain{1}` or `SpatialDomain{1}`.
* `G`: point group. Either `TrivialGroup` or `ReflectionGroup`.
* `W`: wavelet class, e.g. `Morlet` or `Gammatone`.
Its fields are
* `ϕ`: low-pass filter, also called scaling function.
* `ψs`: 3d array of wavelets, indexed by spin, chroma, and octave.
* `behavior`: mutable behavior, e.g. oversampling.
* `spec`: immutable specifications, e.g. number of filters per octave.
To create a `Bank1D`
1. define a `Spec1D`,
2. define a `PathKey` on which to apply the wavelet filterbank,
3. if needed, provide behavior characteristics as keyword arguments.
Example:
spec = Spec1D(nFilters_per_octave = 12, nOctaves = 10)
pathkey = PathKey(:time)
bank = Bank1D(spec, pathkey, j_range = 2:9)
"""
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
            spec::Spec1D{T,D,G,W},
            pathkey::PathKey ;
            is_ϕ_applied::Bool = false,
            j_range::UnitRange{Int} = 0:(spec.nOctaves-1),
            log2_oversampling::Int = 0,
            max_log2_stride::Int = spec.nOctaves-1)
        (nΘs, nΧs, nJs) = size(spec.ψmetas)
        ψs = Array(AbstractFilter{T,D}, (nΘs, nΧs, nJs))
        ψs[1, :, :] =
            pmap(AbstractFilter, spec.ψmetas[1, :, :], fill(spec, nΧs * nJs))
        (nΘs > 1) && spin!(ψs)
        ϕ = AbstractFilter(spec.ϕmeta, spec)
        renormalize!(ϕ, ψs, spec)
        behavior = Behavior(ϕ, ψs, spec,
            is_ϕ_applied, j_range, log2_oversampling, max_log2_stride, pathkey)
        new{T,D,G,W}(ϕ, ψs, behavior, spec)
    end
end
