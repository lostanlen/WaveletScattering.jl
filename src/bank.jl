"""A `Behavior` object contains the mutable information in a filter bank. The
values of these fields may be changed *after* the construction of the filter
bank, without having to recompute the underlying architecture. Fields:

* `γ_range`: range of wavelet indices that are actually used in the
convolutions. Default is all of them (see outer constructor).

* `is_ϕ_applied`: true if and only if the lowpass filter `ϕ` (also known as
scaling function) in addition to the bandpass filters `ψ`'s.

* `log2_oversampling`: base-2 logarithm of the oversampling factor with respect
to the critical sampling rate. Must be positive. Default is 0, i.e. no
oversampling.

* `max_log2_stride`: base-2 logarithm of the maximum distance between
neighboring coefficients after subsampling (also known as *hop size* or
*stride*)"""
type Behavior
    γ_range::UnitRange
    ϕ_log2_sampling::Int
    ψ_log2_samplings::Vector{Int}
    is_ϕ_applied::Bool
    log2_oversampling::Int
    max_log2_stride::Int
end

function Behavior{T}(ϕ::Symmetric1DFilter{T},
        ψs::AbstractArray{AbstractFourier1DFilter{T}},
        spec::AbstractSpec, γ_range::UnitRange, is_ϕ_applied::Bool,
        log2_oversampling::Int, max_log2_stride::Int)
    ϕ_critical_log2_sampling = critical_log2_sampling(ϕ, spec)
    ϕ_log2_sampling =
        max(ϕ_critical_log2_sampling + log2_oversampling, -max_log2_stride)
    ψ_critical_log2_samplings =
        [ critical_log2_sampling(ψ, spec) for ψ in ψs[:,1] ]
    ψ_log2_samplings =
        max(ψ_critical_log2_samplings + log2_oversampling, -max_log2_stride)
    min_log2_sampling = min(ϕ_log2_sampling, minimum(ψ_log2_samplings))
    max_log2_stride = min(max_log2_stride, -min_log2_sampling)
    Behavior(γ_range, ϕ_log2_sampling, ψ_log2_samplings, is_ϕ_applied,
        log2_oversampling, max_log2_stride)
end

"""An `AbstractBank` is a wavelet filter bank. Filter banks are of two kinds:
* `AbstractNonOrientedBank`: no orientation variable `θ`, only a scale
variable `γ`. One-dimensional filter banks with real inputs are non-oriented.
* `AbstractOrientedBank`: wavelets have an orientation variable `θ` and a scale
variable `γ`. Two-dimensional filter banks are oriented. One-dimensional filter
banks with complex outputs are also oriented, in the sense that the sign
of the center frequency is regarded as an orientation variable with two
possible values."""
abstract AbstractBank{T<:Number}
abstract AbstractNonOrientedBank{T<:Number} <: AbstractBank{T}
abstract AbstractOrientedBank{T<:Number} <: AbstractBank{T}

"""A `FourierNonOriented1DBank` is a one-dimensional, non-oriented filter bank
defined in the Fourier domain. It is not oriented in the sense that only the
positive frequencies are guaranteed to be covered. Indeed, assuming that the
input signal will be real — i.e. its Fourier will be symmetric — it can be
recovered from the ""positive half"" of its Fourier spectrum. In summary, it is
advisable to use this type of filter bank when handling real 1d data of
moderate to large length."""
immutable FourierNonOriented1DBank{T<:FFTW.fftwNumber} <:
        AbstractNonOrientedBank{T}
    ϕ::Symmetric1DFilter{T}
    ψs::Vector{AbstractFourier1DFilter{T}}
    behavior::Behavior
    metas::Vector{NonOrientedMeta}
    spec::Abstract1DSpec{T}
    function call{T<:FFTW.fftwNumber}(
            ::Type{FourierNonOriented1DBank{T}}, spec::Abstract1DSpec ;
            γ_range::UnitRange{Int} = 0:typemax(Int),
            is_ϕ_applied::Bool = false, log2_oversampling::Int = 0,
            max_log2_stride::Int = typemax(Int))
        T == spec.signaltype || error("""Type parameter of
        FourierNonOriented1DBankmust must be equal to spec.signaltype""")
        γs, χs, js = gammas(spec), chromas(spec), octaves(spec)
        ξs, qs = centerfrequencies(spec), qualityfactors(spec)
        scs, bws = scales(spec), bandwidths(spec)
        @inbounds metas = [ NonOrientedMeta(
            γs[1+γ], χs[1+γ], bws[1+γ], ξs[1+γ], js[1+γ], qs[1+γ], scs[1+γ])
            for γ in γs ]
        ψs = pmap(fourierwavelet, metas, fill(spec, length(metas)))
        ψs = convert(Array{AbstractFourier1DFilter{T}}, ψs)
        ϕ = scalingfunction(spec)
        renormalize!(ϕ, ψs, metas, spec)
        γ_range = max(γ_range.start, 0):min(γ_range.stop, length(γs) - 1)
        behavior = Behavior(ϕ, ψs, spec,
            γ_range, is_ϕ_applied, log2_oversampling, max_log2_stride)
        new{T}(ϕ, ψs, behavior, metas, spec)
    end
end
FourierNonOriented1DBank(spec::Abstract1DSpec ; args...) =
    FourierNonOriented1DBank{spec.signaltype}(spec ; args)

"""A `FourierOriented1DBank` is a one-dimensional, oriented filter bank defined
in the Fourier domain. It is "oriented" insofar as its filters have negative
center frequencies as well as positive center frequencies, as represented by
the orientation parameter `θ`. It is advisable to use this type of filter bank
when handling complex 1d data of moderate to large length."""
immutable FourierOriented1DBank{T<:FFTW.fftwNumber} <: AbstractOrientedBank{T}
    ϕ::Symmetric1DFilter{T}
    ψs::Matrix{AbstractFourier1DFilter{T}}
    behavior::Behavior
    metas::Matrix{OrientedMeta}
    spec::Abstract1DSpec{T}
    function call{T<:FFTW.fftwNumber}(
            ::Type{FourierOriented1DBank{T}}, spec::Abstract1DSpec ;
            γ_range::UnitRange{Int} = 0:typemax(Int),
            is_ϕ_applied::Bool = false, log2_oversampling::Int = 0,
            max_log2_stride::Int = typemax(Int))
        T == spec.signaltype || error("""Type parameter of
        FourierNonOriented1DBankmust be equal to spec.signaltype""")
        γs, χs, js = gammas(spec), chromas(spec), octaves(spec)
        ξs, qs = centerfrequencies(spec), qualityfactors(spec)
        scs, bws = scales(spec), bandwidths(spec)
        θs = 0:1
        @inbounds metas = [ OrientedMeta(
            γs[γ], θs[θ], χs[γ], bws[γ], ξs[γ], js[γ], qs[γ], scs[γ])
            for γ in eachindex(γs), θ in eachindex(θs) ]
        ψs = pmap(fourierwavelet, metas[:, 1], fill(spec, length(metas)))
        ψs = convert(Array{AbstractFourier1DFilter{T}}, ψs)
        ψs = hcat(ψs, map(spin, ψs))
        ϕ = scalingfunction(spec)
        renormalize!(ϕ, ψs, metas, spec)
        γ_range = max(γ_range.start, 0):min(γ_range.stop, length(γs) - 1)
        behavior = Behavior(ϕ, ψs, spec,
            γ_range, is_ϕ_applied, log2_oversampling, max_log2_stride)
        new{T}(ϕ, ψs, behavior, metas, spec)
    end
end
FourierOriented1DBank(spec::Abstract1DSpec ; args...) =
    FourierOriented1DBank{spec.signaltype}(spec, args...)
