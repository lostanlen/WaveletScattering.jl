abstract AbstractBank{
    T<:Number,
    D<:AbstractDomain,
    G<:AbstractPointGroup,
    W<:RedundantWaveletClass}

immutable NullBank <: AbstractBank
end
