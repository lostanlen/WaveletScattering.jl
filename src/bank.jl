abstract AbstractBank{
    T<:Number,
    D<:AbstractDomain,
    G<:AbstractPointGroup,
    W<:RedundantWaveletClass}

immutable NullBank <: AbstractBank
end

kthrange(syms::Vector{Symbol}, x::AbstractArray, k::Int) =
    (PathKey(syms[k]) => (1:1:size(x, k)))
