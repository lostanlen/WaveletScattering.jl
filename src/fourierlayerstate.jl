# FourierLayerState
immutable FourierLayerState <: AbstractScatteredLayerState
    blobs::Vector{Mocha.Blob}
    blobs_diff::Vector{Mocha.Blob}
    layer::FourierLayer
end

function FourierLayerState(
        backend::Mocha.CPUBackend,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob},
        layer::FourierLayer)
    blobs = Vector{Mocha.Blob}(length(inputs))
    for idblob in eachindex(inputs)
        innodes = inputs[idblob].nodes
        # TODO: tighten typing
        outnodes = DataStructures.SortedDict{Pair,RealFourierNode,
            Base.Order.ForwardOrdering}()
        for path in keys(innodes)
            outnodes[path] = RealFourierNode(
                innodes[path],
                findin(innodes[path].ranges, layer.pathkeys),
                layer.flags,
                layer.timelimit)
        end
        blobs[idblob] = ScatteredBlob(outnodes)
    end
    # TODO: build diffs
    return FourierLayerState(blobs, diffs, layer)
end

function Base.findin{N}(
        ranges::NTuple{N,Pair{PathKey,StepRange{Int,Int}}},
        fourierkeys::Vector{PathKey})
    return find([any(r.first.==fourierkeys) for r in ranges])
end

function Mocha.setup(
        backend::Mocha.Backend,
        layer::FourierLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    return FourierLayerState(backend, diffs, inputs, layer)
end
