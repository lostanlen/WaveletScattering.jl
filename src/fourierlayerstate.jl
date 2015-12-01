# FourierLayerState
immutable FourierLayerState{BANK<:AbstractBank,BLOB<:ScatteredBlob} <:
        AbstractScatteredLayerState
    layer::FourierLayer
    blobs::Vector{BLOB}
    blobs_diff::Any
end
