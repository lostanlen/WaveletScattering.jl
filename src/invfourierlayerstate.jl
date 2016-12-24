# InvFourierLayerState
immutable InvFourierLayerState <: AbstractScatteredLayerState
    blobs::Vector{Mocha.Blob}
    layer::InvFourierLayer
end

function InvFourierLayerState(
        backend::Mocha.CPUBackend,
        layer::InvFourierLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    blobs = Vector{Mocha.Blob}(length(inputs))
    for idblob in eachindex(inputs)
        innodes = inputs[idblob].nodes
        outnodes = DataStructures.SortedDict{Path,InvComplexFourierNode,
            Base.Order.ForwardOrdering}()
        for path in keys(innodes)
            outnodes[path] = InvComplexFourierNode(
                innodes[path],
                findin(innodes[path].ranges, layer.pathkeys),
                layer.flags,
                layer.timelimit)
        end
        blobs[idblob] = ScatteredBlob(outnodes)
    end
    # TODO: build diffs
    return InvFourierLayerState(blobs, layer)
end

function Mocha.setup(
        backend::Mocha.Backend,
        layer::InvFourierLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    return InvFourierLayerState(backend, layer, inputs, diffs)
end
