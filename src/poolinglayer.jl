immutable PoolingLayer <: Mocha.Layer
    bottoms::Vector{Symbol}
    name::AbstractString
    pathkeys::Vector{PathKey}
    pooling::Function
    tops::Vector{Symbol}
end

function PoolingLayer( ;
        bottoms::Vector{Symbol} = Symbol[],
        name::AbstractString = "pointwise",
        pathkeys::Vector{PathKey} = PathKey[],
        pooling::Function = mean,
        tops::Vector{Symbol} = Symbol[])
    @assert length(bottoms) > 0
    @assert length(tops) == length(bottoms)
    PoolingLayer(bottoms, name, pathkeys, pooling, tops)
end

Mocha.can_do_bp(::PoolingLayer) = false
Mocha.has_neuron(::PoolingLayer) = false
Mocha.has_param(::PoolingLayer) = false
