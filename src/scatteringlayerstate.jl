# ScatteredLayerState
immutable ScatteredLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    blobs::Vector{B}
    layer::ScatteringLayer
end
