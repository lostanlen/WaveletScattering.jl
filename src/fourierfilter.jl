abstract AbstractFourierFilter{T<:Number,N} <: AbstractFilter{T,N}
abstract AbstractFourier1DFilter{T<:Number} <: AbstractFourierFilter{T,1}

"""An Analytic1DFilter has only positive frequencies. Its Fourier-domain support is ranges between posfirst and (posfirst+length(pos)-1)<N/2."""
immutable Analytic1DFilter{T<:Number}  <: AbstractFourier1DFilter{T}
    pos::Vector{T}
    posfirst::Int
end

"""A Coanalytic1DFilter has only negative frequencies. Its Fourier-domain support ranges between (neglast-length(neg)+1)>(-N/2) and neglast."""
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

# multiplication operator with scalar *
Base.(:*){T<:Number}(ψ::Analytic1DFilter{T}, b::Number) =
    Analytic1DFilter{T}(b * ψ.pos, ψ.posfirst)
Base.(:*){T<:Number}(ψ::Coanalytic1DFilter{T}, b::Number) =
    Coanalytic1DFilter{T}(b * ψ.neg, ψ.neglast)
Base.(:*){T<:Number}(ψ::FullResolution1DFilter{T}, b::Number) =
    FullResolution1DFilter{T}(b * ψ.coeff)
Base.(:*){T<:Number}(ψ::Vanishing1DFilter{T}, b::Number) =
    Vanishing1DFilter(b * ψ.an, b * ψ.coan)
Base.(:*){T<:Number}(ψ::VanishingWithMidpoint1DFilter{T}, b::Number) =
    VanishingWithMidpoint1DFilter(b * ψ.an, b * ψ.coan, b * ψ.midpoint)
Base.(:*)(b::Number, ψ::AbstractFourierFilter) = ψ * b

# indexing between -N/2 and N/2-1
function Base.getindex{T}(ψ::Analytic1DFilter{T}, i::Integer)
    i<ψ.posfirst && return zero(T)
    i>(ψ.posfirst + length(ψ.pos) - 1) && return zero(T)
    @inbounds return ψ.pos[1 - ψ.posfirst + i]
end
function Base.getindex{T}(ψ::Analytic1DFilter{T}, I::UnitRange{Int64})
    start = max(I.start, ψ.posfirst)
    stop = min(I.stop, ψ.posfirst + length(ψ.pos) - 1)
    return T[
        zeros(T, max(start - I.start, 0));
        ψ.pos[1 - ψ.posfirst + (start:stop)];
        zeros(T, max(I.stop - max(I.start - 1, stop), 0)) ]
end
function Base.getindex{T}(ψ::Coanalytic1DFilter{T}, i::Integer)
    i<(ψ.neglast - length(ψ.neg) + 1) && return zero(T)
    i>ψ.neglast && return zero(T)
    return ψ.neg[i - ψ.neglast + end]
end
function Base.getindex{T}(ψ::Coanalytic1DFilter{T}, I::UnitRange{Int64})
    start = max(I.start, ψ.neglast - length(ψ.neg) + 1)
    stop = min(I.stop, ψ.neglast)
    return T[
        zeros(T, max(min(start, I.stop + 1) - I.start, 0));
        ψ.neg[end - ψ.neglast + (start:stop)];
        zeros(T, max(I.stop - max(I.start - 1, stop), 0)) ]
end
function Base.getindex{T}(ψ::FullResolution1DFilter{T}, i::Integer)
    N = length(ψ.coeff)
    halfN = N >> 1
    i<(-halfN) && return zero(T)
    i>(halfN-1) && return zero(T)
    return ψ.coeff[mod1(1 + i, N)]
end
function Base.getindex{T}(ψ::FullResolution1DFilter{T}, I::UnitRange{Int64})
    N = length(ψ.coeff)
    halfN = N >> 1
    start = max(I.start, -halfN)
    stop = min(I.stop, halfN-1)
    T[
        zeros(T, max(min(start, I.stop + 1) - I.start, 0)) ;
        ψ.coeff[[ mod1(1+ω, N) for ω in start:stop ]] ;
        zeros(T, max(I.stop - max(I.start - 1, stop), 0)) ]
end
function Base.getindex{T}(ψ::Vanishing1DFilter{T}, i::Integer)
    return (i>0 ? ψ.an[i] : ψ.coan[i])
end
function Base.getindex{T}(ψ::Vanishing1DFilter{T}, I::UnitRange{Int64})
    return T[
        ψ.coan[min(0, I.start):min(-1, I.stop)] ;
        zeros(T, Int((I.start < 0) && (I.stop >0))) ;
        ψ.an[max(1, I.start):max(0, I.stop)] ]
end
function Base.getindex(ψ::VanishingWithMidpoint1DFilter, i::Integer)
    halfN = ψ.an.posfirst + length(ψ.an.pos)
    i==(-halfN) && return ψ.midpoint
    return (i>0 ? ψ.an[i] : ψ.coan[i])
end
function Base.getindex{T}(ψ::VanishingWithMidpoint1DFilter{T},
                          I::UnitRange{Int64})
    halfN = ψ.an.posfirst + length(ψ.an.pos)
    output = T[
        ψ.coan[min(0, I.start):min(-1, I.stop)] ;
        ψ.an[max(1, I.start):max(0, I.stop)] ]
    if I.start<=-halfN
        output[-halfN-I.start+1] = ψ.midpoint
    end
    return output
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
function littlewoodpaleyadd!(lp::Vector, ψ::Symmetric1DFilter)
    @inbounds for ω in eachindex(ψ.leg)
        @fastmath temp = abs2(ψ.leg[ω])
        lp[1 + ω] += temp
        lp[N + 1 - ω] += temp
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
wavelet energies (squared magnitudes)."""
function renormalize!{F<:AbstractFourier1DFilter}(ψs::Vector{F},
        ϕ::Symmetric1DFilter, metas::Any, spec::Abstract1DSpec)
    N = 1 << spec.log2_size[1]
    T = spec.signaltype
    if !isinf(spec.max_scale) && spec.max_qualityfactor>1.0
        elbowλ = 1; while (metas[elbowλ].scale<spec.max_scale) elbowλ += 1; end
        elbowω = round(Int, N * metas[elbowλ].centerfrequency)
        λs = elbowλ:length(metas)
        ψmat = zeros(T, (length(ωs), length(λs)))
        for idλ in eachindex(λs) ψmat[:, idλ] = abs2(ψs[λs[idλ]][0:elbowω]); end
        lp = zeros(realtype(T), N)
        for idλ in 1:(elbowλ-1) littlewoodpaleyadd!(lp, ψs[idλ]); end
        littlewoodpaleyadd!(lp, ϕ)
        isa(metas, Vector{NonOrientedMeta}) && symmetrize!(lp)
        remainder = maximum(lp) - lp[1 + (0:elbowω)]
        optimizationmodel = JuMP.Model()
        JuMP.@defVar(optimizationmodel, y[1:length(λs)] >= 0)
        JuMP.@setObjective(optimizationmodel, Min, sum((remainder - ψmat * y)))
        JuMP.@addConstraint(optimizationmodel, remainder .>= ψmat * y)
        JuMP.solve(optimizationmodel)
        ψs[λs] .*= sqrt(2) * JuMP.getValue(y)
    end
    lp = zeros(realtype(T), N)
    for idψ in eachindex(ψs) littlewoodpaleyadd!(lp, ψs[idψ]); end
    isa(metas, Vector{NonOrientedMeta}) && symmetrize!(lp)
    invmax_lp = inv(maximum(lp))
    sqrtinvmax_lp = sqrt(invmax_lp)
    for idψ in eachindex(ψs) ψs[idψ] = ψs[idψ] .* sqrtinvmax_lp; end
    return scale!(lp, invmax_lp)
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
    N = length(lp)
    for ω in 1:(N>>1 - 1)
        halfsum = 0.5 * (lp[1 + ω] + lp[1 + N - ω])
        lp[1 + ω] = halfsum
        lp[1 + N - ω] = halfsum
    end
end
