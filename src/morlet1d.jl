immutable Morlet1DSpec{T<:Number} <: Abstract1DSpec{T}
    ɛ::Float64
    log2_length::Int
    max_qualityfactor::Float64
    max_scale::Float64
    mother_frequency::Float64
    nFilters_per_octave::Int
    nOctaves::Int
    signaltype::Type{T}
    function call{T<:Number}(::Type{Morlet1DSpec{T}}, signaltype::Type{T};
                             ɛ=default_ɛ(T), log2_length=15,
                             max_qualityfactor=nothing, max_scale=Inf,
                             mother_centerfrequency=nothing,
                             nFilters_per_octave=nothing, nOctaves=nothing)
        max_qualityfactor =
            default_max_qualityfactor(max_qualityfactor, nFilters_per_octave)
        nFilters_per_octave =
            default_nFilters_per_octave(max_qualityfactor, nFilters_per_octave)
        motherfrequency =
            default_motherfrequency(Morlet1DSpec, nFilters_per_octave)
        nOctaves =
            default_nOctaves(Morlet1DSpec, nOctaves, log2_length,
                             max_qualityfactor, max_scale, nFilters_per_octave)
        checkspec(ɛ, log2_length, max_qualityfactor,
                  max_scale, nFilters_per_octave, nOctaves)
        new{T}(ɛ, log2_length, max_qualityfactor,
               max_scale, nFilters_per_octave, nOctaves, signaltype)
    end
end

"By default, Morlet1DSpec operates on single-precision real input (Float32)."
Morlet1DSpec(T=Float32; args...) = Morlet1DSpec{T}(T; args...)

function default_nOctaves(::Type{Morlet1DSpec}, nOctaves::Void, log2_length,
                          max_qualityfactor, max_scale, nFilters_per_octave)
    log2_nFilters_per_octave = ceil(Int, log2(nFilters_per_octave))
    log2_max_qualityfactor = ceil(Int, log2(max_qualityfactor))
    if max_scale < (exp2(log2_length)+default_ɛ(T))
        gap = max(1+log2_nFilters_per_octave, 2)
    else
        gap = max(1+log2_nFilters_per_octave, 2+log2_max_qualityfactor)
    end
    nOctaves = log2_length - gap
end

"""
The dimensionless mother center frequency ξ (corresponding to a log-scale γ=0)
is computed as the midpoint between the center frequency of the second center
frequency ξ*2^(-1/nFilters_per_octave) (corresponding to γ=1) and the negative
mother center frequency (1-ξ).

Hence the equation 2ξ = ξ*2^(-1/nFilters_per_octave) + (1-ξ), of which we
derive ξ = 1 / (3.0 - 2^(1/nFilters_per_octave)).

In the special case nFilters_per_octave=1, we manually set ξ=0.39 which is more
accurate with the Littlewood-Paley energy conservation criterion than the
formula ξ=0.4."""
function default_motherfrequency(::Type{Morlet1DSpec}, nFilters_per_octave)
    nFilters_per_octave==1 && return 0.39
    return RealT(inv(3.0 - exp2(-1.0/nFilters_per_octave)))
end

"""Computes in closed form the bandwidths, center frequencies, quality factors,
and scales of all wavelets in a given Morlet 1d filter bank."""
function localize{T<:Number}(spec::Morlet1DSpec{T})
    RealT = realtype(T)
    mother_centerfrequency = spec.nFilters_per_octave==1 ? RealT(0.39) :
        RealT(inv(3.0 - exp2(-1.0/spec.nFilters_per_octave)))
    nΓs = spec.nFilters_per_octave * spec.nOctaves
    resolutions =
        exp2(range(zero(RealT), -one(T)/spec.nFilters_per_octave, nΓs))
    centerfrequencies = mother_centerfrequency * resolutions
    scale_multiplier = sqrt(log(RealT(10.0))/log(RealT(2.0)))
    unbounded_scales =
        scale_multiplier * (spec.max_qualityfactor./centerfrequencies)
    scales = min(unbounded_scales, spec.max_scale)
    unbounded_qualityfactors = scales .* centerfrequencies / scale_multiplier
    qualityfactors = clamp(unbounded_qualityfactors, 1.0, spec.max_qualityfactor)
    bandwidths = resolutions ./ qualityfactors
    scales = scale_multiplier * qualityfactors./centerfrequencies
    return (bandwidths, centerfrequencies, qualityfactors, scales)
end
