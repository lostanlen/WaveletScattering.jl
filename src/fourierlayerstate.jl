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
    nblobs = length(inputs)
    fourierkeys = layer.pathkeys
    for idblob in eachindex(inputs)
        input = inputs[idblob]
        innodes = input.nodes
        outnodes = Dict{Path, AbstractFourierNode}()
        for nodekey in keys(innodes)
            innode = innodes[nodekey]
            region = findin(innode.ranges, fourierkeys)
            outnode = RealFourierNode(
                innodes[nodekey],
                findin(innode.ranges, fourierkeys),
                layer.flags,
                layer.timelimit)
        end
    end
end

function Base.findin{N}(
        ranges::NTuple{N,Pair{PathKey,StepRange{Int,Int}}},
        fourierkeys::Vector{PathKey})
    return find([any(r.first.==fourierkeys) for r in ranges])
end
