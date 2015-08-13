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

function AbstractFourier1DFilter(y, first, last, N)
    if first>(-N/2) && first<0
        if last>0 && last<(N/2)
            # support is in ]-N/2;N/2[
            # we split y between analytic and coanalytic parts
            neg = y[1:(-first)]
            neglast = first
            coan = Coanalytic1DFilter(neg, neglast)
            pos = y[(2-first):end]
            posfirst = 1
            an = Analytic1DFilter(pos, posfirst)
            return Vanishing1DFilter(an, coan)
        elseif last>(N/2-1) && last<N
            # support is in ]-N/2;N[
            # we sum ]-N/2;0[ (or a subset) with ]N/2;N[ (or a subset)
            # we define ]0;N/2[ as the analytic part
            # the midpoint N/2 is handled separately
            neg = vcat(
                y[(2-first+N/2):N],
                y[(1+N):end] + y[1:(1+last-first-N)],
                y[(2+last-first-N):(-first)])
            neglast = -1
            coan = Coanalytic1DFilter(neg, neglast)
            pos = y[(2-first):(-first+N/2)]
            posfirst = 1
            an = Analytic1DFilter(pos, posfirst)
        elseif last>N && last<(3N/2+1)
            # support is in ]-N/2;3N/2]
            # we sum ]-N/2;0[ (or a subset) with ]N/2;N[
            # we sum ]N/2;N[ with ]N;3N/2[ (or a subset)
            # we sum midpoints N/2 and 3N/2 (if present)
            neg = vcat(
                y[(2-first+N/2):N],
                y[(1+N):(N-first)] + y[1:(-first)])
