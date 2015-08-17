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
                y[(2-first):(1+last-first-N)] + y[(2+N-first):(1+last-first)],
                y[(2+last-first-N):(-first+halfN)])
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
The multiplier m is such that:
* `m = 1/max(lp)` for `q = max_q`, i.e. `ξ > max_q*uncertainty/max_s` (ξleft)
* `m = 1/(max_q*max_1)` for `q = 1`, that is `ξ < uncertainty/max_s` (ξright)
* `m` is interpolated linearly in between, that is, for `s=max_s`
If maximum scale is infinite and/or maximum quality factor, the three cases
above collapse into the simpler `m = 1/max(lp)`."""
function renormalize!{T}(ψs, metas, spec::Abstract1DSpec{T})
    N = 1 << spec.log2_size[1]
    lp = zeros(realtype(T), N)
    for λ in eachindex(ψs); littlewoodpaleyadd!(lp, ψs[λ]); end
    if T <: Real
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
        linspaced_qs = linspace(max_qualityfactor, 1, ωright-ωleft+1)
        for ω in 1:(ωleft-1)
            sqrtden = spec.max_qualityfactor
            lp[ω] = lp[ω] / (sqrtden*sqtrden)
        end
        for ω in (ωleft:ωright)
            sqrtden = linspaced_qs[ω-ωleft+1]
            lp[ω] = lp[ω] / (sqrtden*sqrtden)
        end
        invmax_lp = inv(maximum(lp))
        sqrtinvmax_lp = sqrt(invmax_lp)
        centers = round(Int, centerfrequencies(spec)*N)
        for λ in eachindex(ψs)
            ξ = metas[λ].centerfrequency
            ω = 1 + round(Int, N*ξ)
            if ω<ωleft
                normalizer = sqrtinvmax_lp * spec.max_qualityfactor
            elseif ω<ωright
                normalizer = sqrtinvmax_lp * linspaced_qs[ω-ωleft+1]
            else
                normalizer = sqrtinvmax_lp
            end
            ψs[λ] = ψs[λ] / normalizer
        end
    else
        invmax_lp = inv(maximum(lp))
        normalizer = sqrt(invmax_lp)
        for λ in eachindex(ψs); ψs[λ] = ψs[λ] / normalizer; end
    end
    for ω in eachindex(lp); lp[ω] = lp[ω] * invmax_lp; end
    return lp
end

function scalingfunction!{T}(lp::Vector{T}, metas::Vector{AbstractMeta})
    min_ω = round(Int, N * metas[end].centerfrequency)
    phi = [ sqrt(one(T) - lp[1+ω]) for ω in 0:min_ω ]
    sub_last = findlast(phi .< spec.ɛ)
    for ω in 0:sub_last; lp[1+ω] = 1; end
    return Symmetric1DFilter(phi[1+(1:sub_last)], phi[1+0])
end
