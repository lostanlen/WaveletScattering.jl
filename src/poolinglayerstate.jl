immutable PoolingLayerState <: AbstractScatteredLayerState
    blobs::Vector{Mocha.Blob}
    layer::PoolingLayer
end
