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

# element-wise product .*
Base.(:.*){T}(ψ::Analytic1DFilter{T}, x::Number) =
    Analytic1DFilter{T}(ψ.pos .* x, ψ.posfirst)
Base.(:.*){T}(ψ::Analytic1DFilter{T}, x::Vector) =
    Analytic1DFilter{T}(ψ.pos .* x[ψ.posfirst+(1:length(ψ.pos))], ψ.posfirst)
Base.(:.*){T}(ψ::Coanalytic1DFilter{T}, x::Number) =
    Coanalytic1DFilter{T}(ψ.neg .* x, ψ.neglast)
function Base.(:.*){T}(ψ::Coanalytic1DFilter{T}, x::Vector)
    x_range =  ψ.neglast + length(x) + 1 + ((-length(ψ.neg)+1):0)
    Coanalytic1DFilter{T}(ψ.neg .* x[x_range], ψ.neglast)
end
Base.(:.*){T}(ψ::Vanishing1DFilter{T}, x::Union{Number, Vector}) =
    Vanishing1DFilter{T}(ψ.an .* x, ψ.coan .* x)
Base.(:.*){T}(ψ::VanishingWithMidpoint1DFilter{T}, x::Number) =
    VanishingWithMidpoint1DFilter{T}(ψ.an .* x, ψ.coan .* x, midpoint .* x)
Base.(:.*){T}(ψ::VanishingWithMidpoint1DFilter{T}, x::Vector) =
    VanishingWithMidpoint1DFilter{T}(ψ.an .* x, ψ.coan .* x,
                                     ψ.midpoint .* x[1 + end>>1])

# littlewoodpaleyadd!
function littlewoodpaleyadd!(lp::Vector, ψ::Analytic1DFilter)
    @inbounds for ω in eachindex(ψ.an)
        @fastmath lp[ψ.posfirst+ω] += abs2(ψ.an[ω])
    end
end
function littlewoodpaleyadd!(lp::Vector, ψ::Coanalytic1DFilter)
    offset = length(lp) + 1 - length(ψ)
    @inbounds for ω in eachindex(ψ.coan)
        @fastmath lp[offset+ω] += abs2(ψ.coan[ω])
    end
end
function littlewoodpaleyadd!(lp::Vector, ψ::Vanishing1DFilter)
    littlewoodpaleyadd!(lp, ψ.an)
    littlewoodpaleyadd!(lp, ψ.coan)
end
function littlewoodpaleyadd!(lp::Vector, ψ::VanishingWithMidpoint1DFilter)
    littlewoodpaleyadd!(lp, ψ.an)
    littlewoodpaleyadd!(lp, ψ.coan)
    @fastmath lp[end >> 1] += abs2(ψ.midpoint)
end

# littlewoodpaleysum
function littlewoodpaleysum{T}(ψs::Vector{AbstractFourier1DFilter{T}})
    lp = zeros(realtype(T), 1 .<< log2_size[1])
    for λ in eachindex(ψs)
        littlewoodpaleyadd!(lp, ψs[λ])
    end
    for ω in eachindex(lp)
        lp[ω] = sqrt(lp[ω])
    end
end


"""Returns the type parameter of a complex type.
For example, `realtype(Complex{Float32})` returns `Float32`.
For numeric real types, e.g. `Float32`, it is a no-op."""
realtype{T<:Real}(::Type{T}) = T
realtype{T<:Real}(::Type{Complex{T}}) = T

# renormalize!
function renormalize!(ψs, lp, metas, spec::Abstract1DSpec)
    siglength = 1 .<< log2_size[1]
    multiplier = inv(lp)
    for λ in eachindex(ψs)
        ψs[λ] = ψs[λ] .* multiplier
    end
end
