"""A `Behavior` object contains the mutable information in a filter bank. The
values of these fields may be changed *after* the construction of the filter
bank, without having to recompute the underlying architecture. Fields:

* γ_range: range of wavelet indices that are actually used in the convolutions.
Default is all of them (see outer constructor).

* log2_oversampling: base-2 logarithm of the oversampling factor with respect
to the critical sampling rate 2^(-j-1). Must be positive.
Default is 0, i.e. no oversampling.

* min_log2_resolution : base-2 logarithm of the minimal undersampling factor.
Default is (-J+1), i.e. it has no effect. Set min_log2_resolution to 0 to avoid
any undersampling."""
type Behavior
    γ_range::UnitRange
    log2_oversampling::Int
    min_log2_resolution::Int
end
function Behavior(js::Vector)
    γ_range = 0:(length(js)-1)
    log2_oversampling = 0
    min_log2_resolution = -js[end] + 1
    Behavior(γ_range, log2_oversampling, min_log2_resolution)
end

"""An `AbstractBank` is a wavelet filter bank. Filter banks are of two kinds:
* AbstractNonOrientedBank: no orientation variable θ, only a scale variable γ.
One-dimensional filter banks with real inputs are non-oriented.
* AbstractOrientedBank: wavelets have an orientation variable θ and a scale
variable γ. Two-dimensional filter banks are oriented. One-dimensional filter
banks with complex outputs are also oriented, in the sense that the sign
of the center frequency is regarded as an orientation variable with two
possible values."""
abstract AbstractBank{T<:Number}
abstract AbstractNonOrientedBank{T<:Number} <: AbstractBank{T}
abstract AbstractOrientedBank{T<:Number} <: AbstractBank{T}

"""An `FourierNonOriented1DBank` is a one-dimensional, non-oriented filter bank
defined in the Fourier domain."""
immutable FourierNonOriented1DBank{T<:Number} <: AbstractNonOrientedBank{T}
    ψs::Vector{AbstractFourier1DFilter{T}}
    ϕ::Symmetric1DFilter{T}
    behavior::Behavior
    metas::Vector{NonOrientedMeta}
    spec::Abstract1DSpec{T}
    function call{T<:Number}(::Type{FourierNonOriented1DBank{T}},
      spec::Abstract1DSpec)
        T == spec.signaltype || error("""Type parameter of
        FourierNonOriented1DBankmust be equal to spec.signaltype""")
        γs, χs, js = gammas(spec), chromas(spec), octaves(spec)
        ξs, qs = centerfrequencies(spec), qualityfactors(spec)
        scs, bws = scales(spec), bandwidths(spec)
        @inbounds metas = [
            NonOrientedMeta(γs[i], χs[i], bws[i], ξs[i], js[i], qs[i], scs[i])
            for i in eachindex(γs)]
        @inbounds ψs = [fourierwavelet(meta, spec) for meta in metas]
        lp = renormalize!(ψs, metas, spec)
        ϕ = scalingfunction!(lp, metas)
        behavior = Behavior(js)
        new{T}(ψs, ϕ, behavior, metas, spec)
    end
end
FourierNonOriented1DBank(spec::Abstract1DSpec) =
    FourierNonOriented1DBank{spec.signaltype}(spec)

immutable FourierOriented1DBank{T<:Number} <: AbstractOrientedBank{T}
    ψs::Matrix{AbstractFourier1DFilter{T}}
    ϕ::Symmetric1DFilter{T}
    behavior::Behavior
    metas::Matrix{OrientedMeta}
    spec::Abstract1DSpec{T}
    function call{T<:Number}(::Type{FourierOriented1DBank{T}},
      spec::Abstract1DSpec)
        T == spec.signaltype || error("""Type parameter of
        FourierNonOriented1DBankmust be equal to spec.signaltype""")
        γs, χs, js = gammas(spec), chromas(spec), octaves(spec)
        ξs, qs = centerfrequencies(spec), qualityfactors(spec)
        scs, bws = scales(spec), bandwidths(spec)
        θs = 1:2
        @inbounds metas = [
            NonOrientedMeta(γs[i], χs[i], bws[i], ξs[i], js[i], qs[i], scs[i])
            for i in eachindex(γs)]
    end
end
