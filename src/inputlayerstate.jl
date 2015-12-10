immutable InputLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    layer::InputLayer
    blob::B
end
