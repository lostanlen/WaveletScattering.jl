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
            # we define ]0;N/2[ as the analytic part
            # the midpoint N/2 is handled separately
            neg = vcat(
                y[(2-first+halfN):N],
                y[(1+N):end] + y[1:(1+last-first-N)],
                y[(2+last-first-N):(-first)])
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
                y[(1-first):(1+last-first-N)] + y[(1+N-first):(-first+3halfN)],
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
        if last<0 && last<(halfN)
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
            pos = vcat(
                y[(2+N-first):N],
                y[1:(halfN-first)]+y[(1+N):(3halfN-first)],
                y[3halfN-first+1:end])
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
