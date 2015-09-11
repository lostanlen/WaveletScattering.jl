abstract AbstractFourierFilter{T<:Number} <: AbstractFilter{T}
abstract AbstractFourier1DFilter{T<:Number} <: AbstractFourierFilter{T}

immutable Analytic1DFilter{T<:Number}  <: AbstractFourier1DFilter{T}
    pos::Vector{T}
    posfirst::Int
end

immutable Coanalytic1DFilter{T<:Number} <: AbstractFourier1DFilter{T}
    neg::Vector{T}
    neglast::Int
end

immutable FullResolution1DFilter{T<:Number} <: AbstractFourier1DFilter{T}
    coeff::Vector{T}
end

immutable Vanishing1DFilter{T<:Number} <: AbstractFourier1DFilter{T}
    an::Analytic1DFilter
    coan::Coanalytic1DFilter{T}
end

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
Base.(:.*){T<:Number}(ψ::FullResolution1DFilter{T}, b) =
    FullResolution1DFilter{T}(ψ.coeff .* b)
Base.(:.*){T<:Number}(ψ::Analytic1DFilter{T}, b::Number) =
    Analytic1DFilter{T}(ψ.pos .* b, ψ.posfirst)
Base.(:.*){T<:Number}(ψ::Coanalytic1DFilter{T}, b::Number) =
    Coanalytic1DFilter{T}(ψ.neg .* b, ψ.neglast)
Base.(:.*){T<:Number}(ψ::Vanishing1DFilter{T}, b::Number) =
    Vanishing1DFilter(scale(ψ.an, b), scale(ψ.coan, b))
Base.(:.*){T<:Number}(ψ::VanishingWithMidpoint1DFilter{T}, b::Number) =
    VanishingWithMidpoint1DFilter(scale(ψ.an, b), scale(ψ.coan, b), ψ.midpoint*b)
Base.(:.*)(b::Number, ψ::AbstractFourierFilter) = ψ .* b

# right division operator /
Base.(:/){T}(ψ::Analytic1DFilter{T}, x::Number) =
    Analytic1DFilter{T}(ψ.pos / x, ψ.posfirst)
Base.(:/){T}(ψ::Coanalytic1DFilter{T}, x::Number) =
    Coanalytic1DFilter{T}(ψ.neg / x, ψ.neglast)
Base.(:/){T}(ψ::Vanishing1DFilter{T}, x::Number) =
    Vanishing1DFilter{T}(ψ.an / x, ψ.coan / x)
Base.(:/){T}(ψ::VanishingWithMidpoint1DFilter{T}, x::Number) =
    VanishingWithMidpoint1DFilter{T}(ψ.an / x, ψ.coan / x, ψ.midpoint / x)

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
* `b = 1/max(lp)` for `q = max_q`, i.e. `ξ > max_q*uncertainty/max_s` (ξleft)
* `b = 1/(max_q*max_1p)` for `q = 1`, that is `ξ < uncertainty/max_s` (ξright)
* `b` is interpolated linearly in between, that is, for `s=max_s`
If maximum scale is infinite and/or maximum quality factor, the three cases
above collapse into the simpler `m = 1/max(lp)`."""
function renormalize!{F<:AbstractFourier1DFilter}(ψs::Vector{F},
        metas, spec::Abstract1DSpec)
    N = 1 << spec.log2_size[1]
    T = spec.signaltype
    lp = zeros(realtype(T), N)
    for λ in eachindex(ψs); littlewoodpaleyadd!(lp, ψs[λ]); end
    if isa(metas, Vector{NonOrientedMeta})
        for ω in 1:(N>>1 - 1)
            halfsum = 0.5 * (lp[1 + ω] + lp[1 + N - ω])
            lp[1 + ω] = halfsum
            lp[1 + N - ω] = halfsum
        end
    end
    if !isinf(spec.max_scale) && spec.max_qualityfactor>1.0
        ξleft = uncertainty(spec) / spec.max_scale
        ξright = spec.max_qualityfactor * ξleft
        ωleft = 1 + round(Int, N * ξleft)
        ωright = 1 + round(Int, N * ξright)
        inv_max_Q = inv(spec.max_qualtyfactor)
        scale!(lp[1:(ωleft-1)], inv_max_Q)
        linspaced_qs = linspace(spec.max_qualityfactor, 1, ωright-ωleft+1)
        inv_linspaced_qs = 1.0 ./ linspaced_qs
        for ω in (ωleft:ωright)
            lp[ω] *= inv_linspaced_qs[ω-ωleft+1] * inv_linspaced_qs[ω-ωleft+1]
        end
        inv_max_lp = inv(maximum(lp))
        sqrtinv_max_lp = sqrt(inv_max_lp)
        sqrtinv_max_Q = sqrt(inv_max_Q)
        sqrtinv_linspaced_qs = sqrt(inv_linspaced_qs)
        centers = [ 1 + round(Int, meta.centerfrequency*N) for meta in metas ]
        for λ in eachindex(ψs)
            ξ = metas[λ].centerfrequency
            ω = 1 + round(Int, N * ξ)
            if ω<ωleft
                b = sqrtinv_max_lp * sqrtinv_max_Q
            elseif ω<ωright
                b = sqrtinv_max_lp * sqrtinv_linspaced_qs[ω-ωleft+1]
            else
                b = sqrtinv_max_lp
            end
            ψs[λ] = ψs[λ] .* b
        end
    else
        invmax_lp = inv(maximum(lp))
        b = sqrt(invmax_lp)
        for λ in eachindex(ψs); ψs[λ] = ψs[λ] .* b; end
    end
    return scale!(lp, invmax_lp)
end

function scalingfunction!{T, M<:AbstractMeta}(lp::Vector{T}, metas::Vector{M})
    firstpeak = sqrt(metas[end].centerfrequency * metas[end-1].centerfrequency)
    min_ω = round(Int, length(lp) * firstpeak)
    @inbounds phi = [ sqrt(one(T) - lp[1+ω]) for ω in 0:min_ω ]
    for ω in 0:min_ω; @inbounds lp[1+ω] = 1; end
    return Symmetric1DFilter(phi[2:end], phi[1+0])
end

spin(ψ::Analytic1DFilter) = Coanalytic1DFilter(reverse(ψ.posfirst), -ψ.posfirst)
spin(ψ::Coanalytic1DFilter) = Analytic1DFilter(reverse(ψ.neglast), -ψ.neglast)
spin(ψ::FullResolution1DFilter) = FullResolution1DFilter(reverse(ψ.coeff))
spin(ψ::Vanishing1DFilter) = Vanishing1DFilter(reverse(ψ.coan), reverse(ψ.an))
spin(ψ::VanishingWithMidpoint1DFilter) =
    VanishingWithMidpoint1DFilter(reverse(ψ.coan), reverse(ψ.an), midpoint)
