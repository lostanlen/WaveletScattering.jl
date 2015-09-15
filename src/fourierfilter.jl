abstract AbstractFourierFilter{T<:Number,N} <: AbstractFilter{T,N}
abstract AbstractFourier1DFilter{T<:Number} <: AbstractFourierFilter{T,1}

"""An Analytic1DFilter has only positive frequencies. Its Fourier-domain support
is ranges between posfirst and (posfirst+length(pos)-1)<N/2."""
immutable Analytic1DFilter{T<:Number}  <: AbstractFourier1DFilter{T}
    pos::Vector{T}
    posfirst::Int
end

"""A Coanalytic1DFilter has only negative frequencies. Its Fourier-domain
support ranges between (neglast-length(neg)+1)>(-N/2) and neglast."""
immutable Coanalytic1DFilter{T<:Number} <: AbstractFourier1DFilter{T}
    neg::Vector{T}
    neglast::Int
end

"""A FullResolution1DFilter spans the entire spectrum. Its is defined in the
Fourier domain as a vector of size N."""
immutable FullResolution1DFilter{T<:Number} <: AbstractFourier1DFilter{T}
    coeff::Vector{T}
end

"""A Vanishing1DFilter has a limited Fourier-domain support, yet spanning both
in positive and negative frequencies. Moreover, it is null for ω=0 and
at the midpoint ω=N/2. Thus, it is defined as the combination of an
Analytic1DFilter (positive frequencies) and a Coanalytic1DFilter (negative
frequencies)."""
immutable Vanishing1DFilter{T<:Number} <: AbstractFourier1DFilter{T}
    an::Analytic1DFilter
    coan::Coanalytic1DFilter{T}
end

"""A VanishingWithMidpoint1DFilter has a limited Fourier-domain support, yet
spanning both in positive and negative frequencies. It is null for ω=0 but
nonzero at the midpoint ω=N/2. Thus, it is defined as the combination of
an Analytic1DFilter (positive frequencies), a Coanalytic1DFilter (negative
frequencies), and a midpoint (frequency ω=N/2)."""
immutable VanishingWithMidpoint1DFilter{T<:Number} <: AbstractFourier1DFilter{T}
    an::Analytic1DFilter{T}
    coan::Coanalytic1DFilter{T}
    midpoint::T
end

function AbstractFourier1DFilter(y, first, last, log2_length)
    N = 1 << log2_length
    halfN = N >> 1
    if first==(-halfN)
        last==(halfN-1) && return FullResolution1DFilter(fftshift(y))
        midpoint = y[1]
        if last>0
            coan = Coanalytic1DFilter(y[1+(1:(halfN-1))], -1)
            an = Analytic1DFilter(y[(1+halfN+1):end], 1)
        else
            an = Analytic1DFilter(zero(typeof(T)), halfN-1)
            coan = Coanalytic1DFilter(y[(1+1):end], last)
        end
        return VanishingWithMidpoint1DFilter(coan, an, midpoint)
    end
    first>0 && return Analytic1DFilter(y, first)
    last<0 && return Coanalytic1DFilter(y, first)
    an = Analytic1DFilter(y[(1-first):end], 1)
    coan = Coanalytic1DFilter(y[1:(-first)], -1)
    return Vanishing1DFilter(an, coan)
end

# element-wise multiplication operator .*
Base.(:.*){T<:Number}(ψ::Analytic1DFilter{T}, b::Number) =
    Analytic1DFilter{T}(b.*ψ.pos, ψ.posfirst)
Base.(:.*){T<:Number}(ψ::Coanalytic1DFilter{T}, b::Number) =
    Coanalytic1DFilter{T}(b.*ψ.neg, ψ.neglast)
Base.(:.*){T<:Number}(ψ::FullResolution1DFilter{T}, b::Number) =
    FullResolution1DFilter{T}(b.*ψ.coeff)
Base.(:.*){T<:Number}(ψ::Vanishing1DFilter{T}, b::Number) =
    Vanishing1DFilter(b.*ψ.an, b.*ψ.coan)
Base.(:.*){T<:Number}(ψ::VanishingWithMidpoint1DFilter{T}, b::Number) =
    VanishingWithMidpoint1DFilter(b.*ψ.an, b.*ψ.coan, b*ψ.midpoint)
Base.(:.*)(b::Number, ψ::AbstractFourierFilter) = ψ .* b

# indexing between -N/2 and N/2-1
function Base.getindex{T}(ψ::Analytic1DFilter{T}, i::Integer)
    i<ψ.posfirst && return zero(T)
    i>(ψ.posfirst + length(ψ.pos) - 1) && return zero(T)
    @inbounds return ψ.pos[1 - ψ.posfirst + i]
end
function Base.getindex{T}(ψ::Analytic1DFilter{T}, I::UnitRange{Int64})
    start = max(I.start, ψ.posfirst)
    stop = min(I.stop, ψ.posfirst + length(ψ.pos)) - 1
    return [
        zeros(T, max(start-I.start, 0));
        ψ.pos[1 - ψ.posfirst + (start:stop)];
        zeros(T, max(I.stop-stop, 0)) ]
end
function Base.getindex{T}(ψ::Coanalytic1DFilter{T}, i::Integer)
    i<(ψ.neglast - length(ψ.neg) + 1) && return zero(T)
    i>ψ.neglast && return zero(T)
    return ψ.neg[1 - (ψ.neglast - length(ψ.neg)) + i]
end
function Base.getindex{T}(ψ::FullResolution1DFilter{T}, i::Integer)
    halfN = length(ψ.coeff) >> 1
    i<(-halfN) && return zero(T)
    i>(halfN-1) && return zero(T)
    return ψ.coeff[1 + halfN + i]
end
function Base.getindex{T}(ψ::FullResolution1DFilter{T}, I::UnitRange{Int64})
    halfN = length(ψ.coeff) >> 1
    start = max(I.start, -halfN)
    stop = min(I.stop, halfN-1)
    return [
        zeros(T, max(start-I.start, 0));
        ψ.coeff[1 + (start:stop) + halfN];
        zeros(T, max(I.stop-stop, 0)) ]
end
function Base.getindex{T}(ψ::Vanishing1DFilter{T}, i::Integer)
    return (i>0 ? ψ.an[i] : ψ.coan[i])
end
function Base.getindex{T}(ψ::Vanishing1DFilter{T}, I::UnitRange{Int64})
    return [
        ψ.coan[min(0, I.start):min(0, I.stop)];
        ψ.coan[max(0, I.start):max(0, I.start)] ]
end
function Base.getindex(ψ::VanishingWithMidpoint1DFilter, i::Integer)
    halfN = ψ.an.posfirst + length(ψ.an.pos)
    i==halfN && return ψ.midpoint
    return (i>0 ? ψ.an[i] : ψ.coan[length(ψ.coan.neg) - ψ.neglast + i])
end
function Base.getindex(ψ::VanishingWithMidpoint1DFilter, I::UnitRange{Int64})
    halfN = ψ.an.posfirst + length(ψ.pos)
    output = [
        ψ.coan[min(0, I.start):min(0, I.stop)] ;
        ψ.coan[max(0, I.start):max(0, min(I.stop, halfN-1))] ]
    return (I.stop==halfN ? [ψ.midpoint; output] : output)
end

"""Adds the squared magnitude of a Fourier-domain wavelet `ψ` to an accumulator
`lp` (squared Littlewood-Paley sum)."""
function littlewoodpaleyadd!(lp::Vector, ψ::Analytic1DFilter)
    @inbounds for ω in eachindex(ψ.pos)
        @fastmath lp[ψ.posfirst+ω] += abs2(ψ.pos[ω])
    end
end
function littlewoodpaleyadd!(lp::Vector, ψ::Coanalytic1DFilter)
    @inbounds for ω in eachindex(ψ.neg)
        @fastmath lp[length(lp) - length(ψ.neg)+ω] += abs2(ψ.neg[ω])
    end
end
function littlewoodpaleyadd!(lp::Vector, ψ::FullResolution1DFilter)
    @inbounds for ω in eachindex(ψ.coeff)
        @fastmath lp[ω] += abs2(ψ.coeff[ω])
    end
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
Base.maximum(ψ::Analytic1DFilter) = sqrt(maximum(abs2(ψ.pos)))
Base.maximum(ψ::Coanalytic1DFilter) = sqrt(maximum(abs2(ψ.neg)))
Base.maximum(ψ::FullResolution1DFilter) = sqrt(maximum(abs2(ψ.coeff)))
Base.maximum(ψ::Vanishing1DFilter) = max(maximum(ψ.an), maximum(ψ.coan))
Base.maximum(ψ::VanishingWithMidpoint1DFilter) =
    max(maximum(ψ.an), maximum(ψ.coan), abs(ψ.midpoint))

"""Returns the type parameter of a complex type.
For example, `realtype(Complex{Float32})` returns `Float32`.
For numeric real types, e.g. `Float32`, it is a no-op."""
realtype{T<:Real}(::Type{T}) = T
realtype{T<:Real}(::Type{Complex{T}}) = T

"""Renormalizes an array of Fourier-domain wavelets `ψs` by the square root
of the maximum value of its Littlewood-Paley sum `lp`. This ensures approximate
an energy conservation property for the subsequent wavelet transform operator.
The Littlewood-Paley sum `lp` is defined, for each frequency `ω`, as the sum of
wavelet energies (squared magnitudes).
If the quality factor varies across frequencies, the multiplier is no longer a
scalar number, since it adapts to the quality factor q at every frequency ξ.
The multiplier b is such that:
* `b = 1/max(lp)` for `q = max_q`, i.e. `ξ > max_q*uncertainty/max_s` (`ξleft`)
* `b = 1/(max_q*max_1p)` for `q = 1`, that is `ξ < uncertainty/max_s` (`ξright`)
* `b` is interpolated linearly in between, that is, for `s = max_s`
If maximum scale is infinite and/or maximum quality factor, the three cases
above collapse into the simpler `m = 1/max(lp)`."""
function renormalize!{F<:AbstractFourier1DFilter}(ψs::Vector{F},
        metas, spec::Abstract1DSpec)
    N = 1 << spec.log2_size[1]
    T = spec.signaltype
    lp = zeros(realtype(T), N)
    for λ in eachindex(ψs); littlewoodpaleyadd!(lp, ψs[λ]); end
    isa(metas, Vector{NonOrientedMeta}) && symmetrize!(lp)
    if !isinf(spec.max_scale) && spec.max_qualityfactor>1.0
        ξleft = uncertainty(spec) / spec.max_scale
        ξright = spec.max_qualityfactor * ξleft
        λs = find([metas[λ].centerfrequency < ξelbow for λ in eachindex(ψs)])
        ωs = round(Int, N * ξleft):round(Int, N * ξright)
        ψmat = []

        lp = zeros(realtype(T), N)
        for λ in eachindex(ψs); littlewoodpaleyadd!(lp, ψs[λ]); end
        isa(metas, Vector{NonOrientedMeta}) && symmetrize!(lp)
        invmax_lp = inv(max_lp)
    else
        invmax_lp = inv(maximum(lp))
        sqrtinvmax_lp = sqrt(invmax_lp)
        for λ in eachindex(ψs); ψs[λ] = ψs[λ] .* sqrtinvmax_lp; end
    end
    return scale!(lp, invmax_lp)
end

"""Given a Littlewood-Paley sum `lp`, constructs a real symmetric low-pass
filter in the Fourier domain. This low-pass filter is called the ""scaling
function"" of the wavelet filterbank, as it ensures energy conservation
and invertibility."""
function scalingfunction!{T, M<:AbstractMeta}(lp::Vector{T}, metas::Vector{M})
    firstpeak = sqrt(metas[end].centerfrequency * metas[end-1].centerfrequency)
    min_ω = round(Int, length(lp) * firstpeak)
    @inbounds leg = [ sqrt(one(T) - lp[1+ω]) for ω in 1:min_ω ]
    for ω in 0:min_ω; @inbounds lp[1+ω] = 1; end
    return Symmetric1DFilter(leg, one(T))
end

"""Reverses the coefficients of a Fourier-domain 1D filter `ψ` to yield a
""mirror filter"" whose center frequency is of opposite sign. This corresponds
to a conjugation in the time domain."""
spin(ψ::Analytic1DFilter) = Coanalytic1DFilter(reverse(ψ.posfirst), -ψ.posfirst)
spin(ψ::Coanalytic1DFilter) = Analytic1DFilter(reverse(ψ.neglast), -ψ.neglast)
spin(ψ::FullResolution1DFilter) = FullResolution1DFilter(reverse(ψ.coeff))
spin(ψ::Vanishing1DFilter) = Vanishing1DFilter(reverse(ψ.coan), reverse(ψ.an))
spin(ψ::VanishingWithMidpoint1DFilter) =
    VanishingWithMidpoint1DFilter(reverse(ψ.coan), reverse(ψ.an), midpoint)

function symmetrize!(lp::Vector)
    for ω in 1:(length(lp)>>1 - 1)
        halfsum = 0.5 * (lp[1 + ω] + lp[1 + length(lp)-ω])
        lp[1 + ω] = halfsum
        lp[1 + N-ω] = halfsum
    end
end
