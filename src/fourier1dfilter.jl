"""An Analytic1DFilter has only positive frequencies. Its Fourier-domain support
ranges between posfirst and (posfirst+length(pos)-1)<N/2."""
immutable Analytic1DFilter{T} <: AbstractFilter{T,FourierDomain{1}}
    pos::Vector{T}
    posfirst::Int
end

"""A Coanalytic1DFilter has only negative frequencies. Its Fourier-domain
support ranges between (neglast-length(neg)+1)>(-N/2) and neglast."""
immutable Coanalytic1DFilter{T} <: AbstractFilter{T,FourierDomain{1}}
    neg::Vector{T}
    neglast::Int
end

"""A FourierSymmetric1DFilter has a real impulse respnse in the spatial domain,
hence a symmetric response in the Fourier domain. It consists of a leg of
nonzero frequencies and a ""zero"" coefficient."""
immutable FourierSymmetric1DFilter{T} <: AbstractFilter{T,FourierDomain{1}}
    leg::Vector{T}
    zero::T
end

"""A FullResolution1DFilter spans the entire spectrum. Its is defined in the
Fourier domain as a vector of size N."""
immutable FullResolution1DFilter{T} <: AbstractFilter{T,FourierDomain{1}}
    coeff::Vector{T}
end

"""A Vanishing1DFilter has a limited Fourier-domain support, yet spanning both
in positive and negative frequencies. Moreover, it is null for ω=0 and
at the midpoint ω=N/2. Thus, it is defined as the combination of an
Analytic1DFilter (positive frequencies) and a Coanalytic1DFilter (negative
frequencies)."""
immutable Vanishing1DFilter{T} <: AbstractFilter{T,FourierDomain{1}}
    an::Analytic1DFilter
    coan::Coanalytic1DFilter{T}
end

"""A VanishingWithMidpoint1DFilter has a limited Fourier-domain support, yet
spanning both in positive and negative frequencies. It is null for ω=0 but
nonzero at the midpoint ω=N/2. Thus, it is defined as the combination of
an Analytic1DFilter (positive frequencies), a Coanalytic1DFilter (negative
frequencies), and a midpoint (frequency ω=N/2)."""
immutable VanishingWithMidpoint1DFilter{T} <: AbstractFilter{T,FourierDomain{1}}
    an::Analytic1DFilter{T}
    coan::Coanalytic1DFilter{T}
    midpoint::T
end

"""Given an array `y` and a filter bank specification `spec`, the
constructor `AbstractFilter` truncates `y` according to the magnitude
threshold `spec.ɛ` and returns the appropriate subtype of `AbstractFilter`
to represent the truncated array. The output is cast with the supertype
`AbstractFilter` so that the output of this function can be stored
efficiently into an array of abstract objects.
For one-dimensional filters in the Fourier domain, the output type may be:
* `Analytic1DFilter` if `y` is non-negligible only for positive frequencies
* `Coanalytic1DFilter` if `y` is non-negligible only for negative frequencies
* `Vanishing1DFilter` if the support of `y` encompasses a portion of
  positive frequencies and a portion of negative frequencies, excluding
  midpoint.
* `VanishingWithMidpoint1DFilter` if the support of `y` encompasses a portion
  of positive frequencies and a portion of negative frequencies, including
  midpoint.
* `FullResolution1DFilter` if `y` is non-negligible along the whole spectrum"""
function AbstractFilter{T}(y::Vector{T}, spec::AbstractSpec{T,FourierDomain{1}})
    supertype = AbstractFilter{T,FourierDomain{1}}
    N = 1 << spec.log2_size
    halfN = N >> 1
    ɛ2 = T(spec.ɛ * spec.ɛ)
    y2 = abs2.(y)
    negbools = y2[2:halfN] .> ɛ2
    negfirst, neglast = findfirst(negbools), findlast(negbools)
    posbools = y2[(2+halfN):end] .> ɛ2
    posfirst, poslast = findfirst(posbools), findlast(posbools)
    hasmidpoint = y2[1] .> ɛ2
    if hasmidpoint
        (negfirst == 1) && (neglast == (halfN-1)) &&
            (posfirst == 1) && (poslast == (halfN-1)) &&
            return FullResolution1DFilter(fftshift(y))::supertype
        midpoint = y[1]
        if (neglast == 0)
            coan = Coanalytic1DFilter(zeros(T, 1), -halfN + 1)
        else
            coan = Coanalytic1DFilter(y[(1+negfirst):(1+neglast)],
                neglast - halfN)
        end
        if (poslast == 0)
            an = Analytic1DFilter(zeros(T, 1), halfN - 1)
        else
            an = Analytic1DFilter(y[(1+halfN+posfirst):(1+halfN+poslast)],
                posfirst)
        end
        return VanishingWithMidpoint1DFilter(an, coan, midpoint)::supertype
    end
    if (negfirst == 0)
        return Analytic1DFilter(y[(1+halfN+posfirst):(1+halfN+poslast)],
            posfirst)::supertype
    end
    if (posfirst == 0)
        return Coanalytic1DFilter(y[(1+negfirst):(1+neglast)],
            neglast)::supertype
    end
    an = Analytic1DFilter(y[(1+halfN+posfirst):(1+halfN+poslast)], posfirst)
    coan = Coanalytic1DFilter(y[(1+negfirst):(1+neglast)], neglast - halfN)
    return Vanishing1DFilter(an, coan)::supertype
end

# multiplication operator with scalar *
Base.:*{T}(ψ::Analytic1DFilter{T}, b::Number) =
    Analytic1DFilter{T}(b * ψ.pos, ψ.posfirst)
Base.:*{T}(ψ::Coanalytic1DFilter{T}, b::Number) =
    Coanalytic1DFilter{T}(b * ψ.neg, ψ.neglast)
Base.:*{T}(ψ::FullResolution1DFilter{T}, b::Number) =
    FullResolution1DFilter{T}(b * ψ.coeff)
Base.:*{T}(ψ::FourierSymmetric1DFilter{T}, b::Number) =
    FourierSymmetric1DFilter{T}(b * ψ.leg, T(b) * ψ.zero)
Base.:*{T}(ψ::Vanishing1DFilter{T}, b::Number) =
    Vanishing1DFilter{T}(b * ψ.an, b * ψ.coan)
Base.:*{T}(ψ::VanishingWithMidpoint1DFilter{T}, b::Number) =
    VanishingWithMidpoint1DFilter{T}(b * ψ.an, b * ψ.coan, T(b) * ψ.midpoint)

"""Given a narrowband filter `ψ` and its corresponding filter bank
specification `spec`, returns the base-2 logarithm of the resampling factor
such that the passband of `ψ` is not aliased. This resmpling factor is
expressed in fraction of signal length ; thus, its base-2 logarithm is
non-positive. The result does not account for manual oversampling."""
critical_log2_sampling(ψ::Analytic1DFilter, spec::AbstractSpec) =
    1 + nextpow2_exponent(ψ.posfirst + length(ψ.pos) - 1) - spec.log2_size
critical_log2_sampling(ψ::Coanalytic1DFilter, spec::AbstractSpec) =
    1 + nextpow2_exponent(ψ.neglast - length(ψ.neg) + 1) - spec.log2_size
critical_log2_sampling(ψ::FullResolution1DFilter, spec::AbstractSpec) = 0
critical_log2_sampling(ψ::FourierSymmetric1DFilter, spec::AbstractSpec) =
    1 + nextpow2_exponent(length(ψ.leg)) - spec.log2_size
critical_log2_sampling(ψ::Vanishing1DFilter, spec::AbstractSpec) = max(
    critical_log2_sampling(ψ.an, spec), critical_log2_sampling(ψ.coan, spec))
critical_log2_sampling(ψ::VanishingWithMidpoint1DFilter, spec::AbstractSpec) = 0

# indexing between -N/2 and N/2-1
function Base.getindex{T}(ψ::Analytic1DFilter{T}, i::Integer)
    i < ψ.posfirst && return zero(T)
    i > (ψ.posfirst + length(ψ.pos) - 1) && return zero(T)
    @inbounds return ψ.pos[1 - ψ.posfirst + i]
end
function Base.getindex{T}(ψ::Analytic1DFilter{T}, I::UnitRange{Int64})
    I.stop < ψ.posfirst && return zeros(T, length(I))
    start = max(I.start, ψ.posfirst)
    stop = min(I.stop, ψ.posfirst + length(ψ.pos) - 1)
    return T[
        zeros(T, max(start - I.start, 0)) ;
        ψ.pos[1 - ψ.posfirst + (start:stop)] ;
        zeros(T, max(I.stop - max(I.start - 1, stop), 0)) ]
end
function Base.getindex{T}(ψ::Coanalytic1DFilter{T}, i::Integer)
    i<(ψ.neglast - length(ψ.neg) + 1) && return zero(T)
    i>ψ.neglast && return zero(T)
    return ψ.neg[i - ψ.neglast + end]
end
function Base.getindex{T}(ψ::Coanalytic1DFilter{T}, I::UnitRange{Int64})
    I.start > ψ.neglast && return zeros(T, length(I))
    start = max(I.start, ψ.neglast - length(ψ.neg) + 1)
    stop = min(I.stop, ψ.neglast)
    return T[
        zeros(T, max(min(start, I.stop + 1) - I.start, 0));
        ψ.neg[end - ψ.neglast + (start:stop)];
        zeros(T, max(I.stop - max(I.start - 1, stop), 0)) ]
end
function Base.getindex{T}(ϕ::FourierSymmetric1DFilter{T}, i::Integer)
    i==0 && return ϕ.zero
    abs(i)>length(ϕ.leg) && return zero(T)
    return ϕ.leg[abs(i)]
end
function Base.getindex{T}(ϕ::FourierSymmetric1DFilter{T}, I::UnitRange{Int64})
    T[ ϕ[ω] for ω in I ]
end
function Base.getindex{T}(ψ::FullResolution1DFilter{T}, i::Integer)
    N = length(ψ.coeff)
    halfN = N >> 1
    if i < -halfN
        return zero(T)
    elseif i < 0
        return ψ.coeff[1+i+N]
    elseif i < halfN
        return ψ.coeff[1+i]
    else
        return zero(T)
    end
end
function Base.getindex{T}(ψ::FullResolution1DFilter{T}, I::UnitRange{Int64})
    N = length(ψ.coeff)
    halfN = N >> 1
    I.start > halfN && return zeros(T, length(I))
    I.stop < -halfN && return zeros(T, length(I))
    return T[
        zeros(T, max(-halfN-I.start, 0));
        ψ.coeff[1 + N + (max(I.start, -halfN):min(I.stop, -1))];
        ψ.coeff[1 + (max(I.start, 0):min(I.stop, halfN-1))];
        zeros(T, max(I.stop-halfN+1, 0))
    ]
end
function Base.getindex{T}(ψ::Vanishing1DFilter{T}, i::Integer)
    return (i > 0 ? ψ.an[i] : ψ.coan[i])
end
function Base.getindex{T}(ψ::Vanishing1DFilter{T}, I::UnitRange{Int64})
    return T[
        ψ.coan[min(0, I.start):min(-1, I.stop)] ;
        zeros(T, Int((I.start < 0) && (I.stop > 0))) ;
        ψ.an[max(1, I.start):max(0, I.stop)] ]
end
function Base.getindex(ψ::VanishingWithMidpoint1DFilter, i::Integer)
    halfN = ψ.an.posfirst + length(ψ.an.pos)
    i==(-halfN) && return ψ.midpoint
    return (i > 0 ? ψ.an[i] : ψ.coan[i])
end
function Base.getindex{T}(ψ::VanishingWithMidpoint1DFilter{T},
                          I::UnitRange{Int64})
    halfN = ψ.an.posfirst + length(ψ.an.pos)
    output = T[
        ψ.coan[min(0, I.start):min(-1, I.stop)] ;
        zeros(T, Int((I.start < 0) && (I.stop > 0))) ;
        ψ.an[max(1, I.start):max(0, I.stop)] ]
    if I.start <= -halfN
        output[-halfN-I.start+1] = ψ.midpoint
    end
    return output
end

"""Adds the squared magnitude of a Fourier-domain wavelet `ψ` to an accumulator
`lp` (squared Littlewood-Paley sum)."""
function littlewoodpaleyadd!(lp::Vector, ψ::Analytic1DFilter)
    @inbounds for ω in eachindex(ψ.pos)
        @fastmath lp[ψ.posfirst + ω] += abs2(ψ.pos[ω])
    end
end
function littlewoodpaleyadd!(lp::Vector, ψ::Coanalytic1DFilter)
    @inbounds for ω in eachindex(ψ.neg)
        @fastmath lp[1 + length(lp) + ψ.neglast - length(ψ.neg) + ω] +=
            abs2(ψ.neg[ω])
    end
end
function littlewoodpaleyadd!(lp::Vector, ψ::FullResolution1DFilter)
    @inbounds for ω in eachindex(ψ.coeff)
        @fastmath lp[ω] += abs2(ψ.coeff[ω])
    end
end
function littlewoodpaleyadd!(lp::Vector, ψ::FourierSymmetric1DFilter)
    @inbounds for ω in eachindex(ψ.leg)
        @fastmath lp[1 + ω] += abs2(ψ.leg[ω])
        @fastmath lp[length(lp) + 1 - ω] += abs2(ψ.leg[ω])
    end
    @fastmath lp[1 + 0] += abs2(ψ.zero)
end
function littlewoodpaleyadd!(lp::Vector, ψ::Vanishing1DFilter)
    littlewoodpaleyadd!(lp, ψ.an)
    littlewoodpaleyadd!(lp, ψ.coan)
end
function littlewoodpaleyadd!(lp::Vector, ψ::VanishingWithMidpoint1DFilter)
    littlewoodpaleyadd!(lp, ψ.an)
    littlewoodpaleyadd!(lp, ψ.coan)
    @fastmath lp[1 + (length(lp)>>1)] += abs2(ψ.midpoint)
end

"""Returns the maximum Fourier-domain absolute value of a filter."""
Base.maximum(ψ::Analytic1DFilter) = sqrt(maximum(abs2.(ψ.pos)))
Base.maximum(ψ::Coanalytic1DFilter) = sqrt(maximum(abs2.(ψ.neg)))
Base.maximum(ψ::FullResolution1DFilter) = sqrt(maximum(abs2.(ψ.coeff)))
Base.maximum(ψ::FourierSymmetric1DFilter) =
    sqrt(max(maximum(abs2.(ψ.leg)), abs2(ψ.zero)))
Base.maximum(ψ::Vanishing1DFilter) = max(maximum(ψ.an), maximum(ψ.coan))
Base.maximum(ψ::VanishingWithMidpoint1DFilter) =
    max(maximum(ψ.an), maximum(ψ.coan), abs(ψ.midpoint))

"""The exponent of the smallest power of two not less than input integer `n`.
Returns `0` for `n==0`, and returns `-nextpow2_exponent(-n)` for negative
arguments."""
nextpow2_exponent(n::Unsigned) = (sizeof(n)<<3)-leading_zeros(n-1)
function nextpow2_exponent(n::Integer)
    if n == 0
        return oftype(n, 0)
    elseif n < 0
        return oftype(n, -nextpow2_exponent(unsigned(-n)))
    else
        return oftype(n, nextpow2_exponent(unsigned(n)))
    end
end

"""Renormalizes an array of Fourier-domain wavelets `ψs` by the square root
of the maximum value of its Littlewood-Paley sum `lp`. This ensures approximate
an energy conservation property for the subsequent wavelet transform operator.
The Littlewood-Paley sum `lp` is defined, for each frequency `ω`, as the sum of
wavelet energies (squared magnitudes)."""
function renormalize!{T<:Number,G<:LineGroups}(
        ϕ::FourierSymmetric1DFilter{T},
        ψs::Array{AbstractFilter{T,FourierDomain{1}},3},
        spec::AbstractSpec{T,FourierDomain{1},G})
    lp = zeros(real(T), 1<<spec.log2_size)
    littlewoodpaleyadd!(lp, ϕ)
    for idψ in eachindex(ψs[1, :, :]) littlewoodpaleyadd!(lp, ψs[idψ]); end
    get_n_orientations(spec.pointgroup)==1 && symmetrize!(lp)
    ψs .*= inv(sqrt(lp[1+end>>1]))
    return scale!(lp, inv(lp[1+end>>1]))
end

"""Reverses the coefficients of a Fourier-domain 1D filter `ψ` to yield a
""mirror filter"" whose center frequency is of opposite sign. This corresponds
to a conjugation in the time domain."""
spin(ψ::Analytic1DFilter) = Coanalytic1DFilter(reverse(ψ.pos), -ψ.posfirst)
spin(ψ::Coanalytic1DFilter) = Analytic1DFilter(reverse(ψ.neg), -ψ.neglast)
spin(ψ::FullResolution1DFilter) = FullResolution1DFilter(reverse(ψ.coeff))
spin(ψ::Vanishing1DFilter) = Vanishing1DFilter(spin(ψ.coan), spin(ψ.an))
spin(ψ::VanishingWithMidpoint1DFilter) =
    VanishingWithMidpoint1DFilter(spin(ψ.coan), spin(ψ.an), ψ.midpoint)

"""Returns the symmetric part of a input vector. This function is called by
`renormalize!` to symmetrize the Littlewood-Paley sum of a one-dimensional
filter bank."""
function symmetrize!(lp::Vector)
    N = length(lp)
    for ω in 1:(N>>1 - 1)
        halfsum = 0.5 * (lp[1 + ω] + lp[1 + N - ω])
        lp[1 + ω] = halfsum
        lp[1 + N - ω] = halfsum
    end
end
