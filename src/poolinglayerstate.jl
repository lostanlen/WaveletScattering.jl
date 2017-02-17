immutable PoolingLayerState <: AbstractScatteredLayerState
    blobs::Vector{Mocha.Blob}
    layer::PoolingLayer
end

function Mocha.setup(
        backend::Mocha.Backend,
        layer::PoolingLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    blobs = Vector{ScatteredBlob}(length(inputs))
    for idblob in eachindex(inputs)
        innodes = inputs[idblob].nodes
        outnodes = DataStructures.SortedDict{Path,AbstractNode,
            Base.Order.ForwardOrdering}()
        for inpath in keys(innodes)
            innode = innodes[inpath]
            inranges = innode.ranges
            inpathkeys = [ pair.first for pair in collect(inranges) ]
            subscripts = vcat(
                [find(pathkey.==inpathkeys) for pathkey in layer.pathkeys]...)
            indata = innode.data
            outdata = layer.pooling(indata, subscripts)
            outranges = collect(inranges)
            for subscript in subscripts
                outranges[subscript] =
                    Pair(inranges[subscript].first,
                        StepRange(
                        inranges[subscript].second.start,
                        inranges[subscript].second.stop -
                        inranges[subscript].second.start + 1,
                        inranges[subscript].second.stop))
            end
            outranges = tuple(outranges...)
            outnode = Node(outdata, outranges)
            outpath = inpath
            outnodes[outpath] = outnode
        end
        blobs[idblob] = ScatteredBlob(outnodes)
    end
    PoolingLayerState(blobs, layer)
end
