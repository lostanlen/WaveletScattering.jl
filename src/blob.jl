# ScatteredBlob
immutable ScatteredBlob{T<:Number, N} <: Mocha.Blob{T, N}
    nodes::Dict{Path, AbstractNode{T, N}}
    subscripts::NTuple{PathKey}
end
