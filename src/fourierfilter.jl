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
    if first>(-halfN) && first<0
        if last>0 && last<(halfN)
            # support is in ]-N/2;N/2[
            # we split y between analytic and coanalytic parts
            neg = y[1:(-first)]
            neglast = first
            coan = Coanalytic1DFilter(neg, neglast)
            pos = y[(2-first):end]
            posfirst = 1
            an = Analytic1DFilter(pos, posfirst)
            return Vanishing1DFilter(an, coan)
        elseif last>(halfN-1) && last<N
            # support is in ]-N/2;N[
            # we sum ]-N/2;0[ (or a subset) with ]N/2;N[ (or a subset)
            # if the subsets are not overlapping, we zero-pad in between
            # we define ]0;N/2[ as the analytic part
            # the midpoint N/2 is handled separately
            neg = vcat(
                y[(2-first+halfN):min(N,length(y))],
                zeros(eltype(y), N-min(N,length(y))),
                y[1+min(N,length(y)):end] + y[1:(1+last-first-N)],
                y[max(1,2+last-first-N):(-first)])
            neglast = -1
            coan = Coanalytic1DFilter(neg, neglast)
            pos = y[(2-first):(-first+halfN)]
            posfirst = 1
            an = Analytic1DFilter(pos, posfirst)
            midpoint = y[1-first+halfN]
            return VanishingWithMidpoint1DFilter(an, coan, midpoint)
        elseif last>N && last<(3halfN+1)
            # support is in ]-N/2;3N/2]
            # we sum ]-N/2;0[ (or a subset) with ]N/2;N[
            # we sum ]N/2;N[ with ]N;3N/2[ (or a subset)
            # we sum midpoints N/2 and 3N/2 (if present)
            neg = vcat(
                y[(2-first+halfN):N],
                y[(1+N):(N-first)] + y[1:(-first)])
            neglast = -1
            coan = Coanalytic1DFilter(neg, neglast)
            pos = vcat(
                y[(2-first):(last-first-N)] + y[(2+N-first):(last-first)],
                y[(2+last-first-N):(-first+halfN-1)])
            posfirst = 1
            an = Analytic1DFilter(pos, posfirst)
            if last == 3halfN
                midpoint = y[1-first+halfN] + y[1-first+3halfN]
            else
                midpoint = y[1-first+halfN]
            end
            return VanishingWithMidpoint1DFilter(an, coan, midpoint)
        end
    elseif first>0 && first<(halfN)
        if last>0 && last<(halfN)
            # support is in ]0;N/2[
            # we just define y as the analytic part
            pos = y
            posfirst = first
            return Analytic1DFilter(pos, posfirst)
        elseif last>(halfN+1) && last<N
            # support is in ]0;N[
            # we split y between analytic and coanalytic parts
            # the midpoint halfN is handled separately
            pos = y[1:(halfN-first)]
            posfirst = first
            an = Analytic1DFilter(pos, posfirst)
            neg = y[(2+halfN-first):end]
            neglast = last - N
            coan = Coanalytic1DFilter(neg, neglast)
            midpoint = y[1+halfN-first]
            return VanishingWithMidpoint1DFilter(an, coan, midpoint)
        elseif last>N && last<(3halfN+1)
            # support is in ]0;3N/2]
            # we sum ]0;N/2[ (or a subset) with ]N;N/2[ (or a subset)
            # we define ]N/2;N[ as the coanalytic part
            # the midpoint N/2 is handled separately
            m = min(N, length(y))
            pos = vcat(
                y[(2+N-first):m],
                y[(1+m-N):(1-first+last-N)] + y[(1+m):(1-first+last)],
                zeros(Int, max(first+N-last-1,0)),
                y[max(1,2-first+last-N):(halfN-first)])
            posfirst = 1
            an = Analytic1DFilter(pos, posfirst)
            neg = y[(2+halfN-first):(N-first)]
            neglast = -1
            coan = Coanalytic1DFilter(neg, neglast)
            if last == 3halfN
                midpoint = y[1+halfN-first] + y[1+3halfN-first]
            else
                midpoint = y[1+halfN-first]
            end
            return VanishingWithMidpoint1DFilter(an, coan, midpoint)
        end
    end
end

# scale
Base.scale{T<:Number}(ψ::Analytic1DFilter{T}, b::Number) =
    Analytic1DFilter{T}(ψ.pos .* b, ψ.posfirst)
Base.scale{T<:Number}(ψ::Coanalytic1DFilter{T}, b::Number) =
    Coanalytic1DFilter{T}(ψ.neg .* b, ψ.neglast)
Base.scale{T<:Number}(ψ::Vanishing1DFilter{T}, b::Number) =
    Vanishing1DFilter(scale(ψ.an, b), scale(ψ.coan, b))
Base.scale{T<:Number}(ψ::VanishingWithMidpoint1DFilter{T}, b::Number) =
    VanishingWithMidpoint1DFilter(scale(ψ.an, b), scale(ψ.coan, b), ψ.midpoint*b)

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
function littlewoodpaleyadd!(lp::Vector, ψ::Vanishing1DFilter)
    littlewoodpaleyadd!(lp, ψ.an)
    littlewoodpaleyadd!(lp, ψ.coan)
end
function littlewoodpaleyadd!(lp::Vector, ψ::VanishingWithMidpoint1DFilter)
    littlewoodpaleyadd!(lp, ψ.an)
    littlewoodpaleyadd!(lp, ψ.coan)
    @fastmath lp[1 + (length(lp) >> 1)] += abs2(ψ.midpoint)
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
        for ω in 1:(N>>1-1)
            halfsum = 0.5 * (lp[1 + ω] + lp[1 + N - ω])
            lp[1+ω] = halfsum
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
            ψs[λ] = scale(ψs[λ], b)
        end
    else
        invmax_lp = inv(maximum(lp))
        b = sqrt(invmax_lp)
        for λ in eachindex(ψs); ψs[λ] = scale(ψs[λ], b); end
    end
    return scale!(lp, invmax_lp)
end

function scalingfunction!{T, M<:AbstractMeta}(lp::Vector{T}, metas::Vector{M})
    firstpeak = sqrt(metas[end].centerfrequency * metas[end-1].centerfrequency)
    min_ω = round(Int, length(lp) * firstpeak)
    phi = [ sqrt(one(T) - lp[1+ω]) for ω in 0:min_ω ]
    for ω in 0:min_ω; lp[1+ω] = 1; end
    return Symmetric1DFilter(phi[2:end], phi[1+0])
end

spin(ψ::Analytic1DFilter) = Coanalytic1DFilter(reverse(ψ.posfirst), -ψ.posfirst)
spin(ψ::Coanalytic1DFilter) = Analytic1DFilter(reverse(ψ.neglast), -ψ.neglast)
spin(ψ::Vanishing1DFilter) = Vanishing1DFilter(reverse(ψ.coan), reverse(ψ.an))
spin(ψ::VanishingWithMidpoint1DFilter) =
    VanishingWithMidpoint1DFilter(reverse(ψ.coan), reverse(ψ.an), midpoint)
