# FourierLayerState
immutable FourierLayerState{B<:ScatteredBlob} <:
        AbstractScatteredLayerState
    layer::FourierLayer
    blobs::Vector{B}
    blobs_diff::Vector{B}
end
