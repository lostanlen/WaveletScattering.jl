# FourierLayerState
immutable FourierLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    layer::FourierLayer
    blobs::Vector{B}
    blobs_diff::Vector{B}
end

function FourierLayerState{B<:ScatteredBlob}(
        backend::Mocha.CPUBackend,
        layer::FourierLayer,
        inputs::Vector{B})
    blobs = Vector{B}(length(inputs))
    for idblob in eachindex(inputs)
        innodes = inputs[idblob].nodes
        outnodes = Dict{Path, AbstractFourierNode}()
        for path in keys(innodes)
            outnodes[path] = RealFourierNode(
                innodes[path],
                findin(innodes[path].ranges, layer.pathkeys),
                layer.flags,
                layer.timelimit)
        end
        blobs[idblob] = ScatteredBlob(outnodes)
    end
    return FourierLayerState(layer, blobs, ScatteredBlob[])
end

function Base.findin{N}(
        ranges::NTuple{N,Pair{PathKey,StepRange{Int,Int}}},
        fourierkeys::Vector{PathKey})
    return find([any(r.first.==fourierkeys) for r in ranges])
end
