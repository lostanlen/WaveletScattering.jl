# FourierLayerState
immutable FourierLayerState <: AbstractScatteredLayerState
    blobs::Vector{Mocha.Blob}
    blobs_diff::Vector{Mocha.Blob}
    layer::FourierLayer
end

function FourierLayerState(
        backend::Mocha.CPUBackend,
        diffs::Vector{Mocha.Blob},
        inputs::Vector{Mocha.Blob},
        layer::FourierLayer)
    blobs = Vector{Mocha.Blob}(length(inputs))
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
        diffs::Vector{Mocha.Blob},
        inputs::Vector{Mocha.Blob},
        layer::FourierLayer)
    return FourierLayerState(backend, diffs, inputs, layer)
end
