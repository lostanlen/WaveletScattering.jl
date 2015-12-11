# FourierLayerState
immutable FourierLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    layer::FourierLayer
    blobs::Vector{B}
    blobs_diff::Vector{B}
end

function FourierLayerState(
        backend::Backend,
        layer::FourierLayer,
        inputs::Vector{Blob})
    nblobs = length(inputs)
    fourierkeys = layer.pathkeys
    #for idblob in eachindex(inputs)
        input = inputs[idblob]
        innodes = input.nodes
        outnodes = Dict{Path, AbstractFourierNode}()
        for nodekey in keys(nodes)
            innode = innodes[nodekey]
            region = findin(ranges, fourierkeys)
            outnode = RealFourierNode(
                innodes[nodekey],
                findin(ranges, fourierkeys),
                layer.flags,
                layer.timelimit)
        end
    #end
end

function Base.findin{N}(
        ranges::NTuple{N,Pair{PathKey,StepRange{Int,Int}}},
        fourierkeys::Vector{PathKey})
    return find([any(r.first.==fourierkeys) for r in ranges])
end
