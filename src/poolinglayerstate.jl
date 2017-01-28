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
            subscripts = findin(inpathkeys, layer.pathkeys)
            indata = innode.data
            outdata = layer.pooling(indata, subscripts)
            outranges = inranges
            for subscript in subscripts
                outranges[subscript].second.step =
                    inranges[subscript].second.stop -
                    inranges[subscript].second.start + 1
            end
            outnode = Node(outdata, outranges)
            outpath = inpath
            outnodes[outpath] = outnode
        end
        blobs[idblob] = ScatteredBlob(outnodes)
    end
    WaveletLayerState(blobs, layer)
end
