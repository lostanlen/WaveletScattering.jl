# WaveletLayerState
immutable WaveletLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    blobs::Vector{B}
    layer::WaveletLayer
end

function Mocha.setup(
        backend::Mocha.Backend,
        layer::WaveletLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    pathkey = layer.bank.behavior.pathkey
    chromakey = prepend(:χ, layer.bank.behavior.pathkey)
    chromarange = chromakey => 0:1:(layer.bank.spec.n_filters_per_octave-1)
    octavekey = prepend(:j, layer.bank.behavior.pathkey)
    innodes = inputs[1].nodes
    outnodes = DataStructures.SortedDict{Path,AbstractNode,
        Base.Order.ForwardOrdering}()
    for inpath in keys(innodes)
        innode = innodes[inpath]
        inranges = innode.ranges
        inpathkeys = [ pair.first for pair in collect(inranges) ]
        subscripts = find(pathkey .== inpathkeys)
        insizes = collect(size(innode.data))
        if isa(innode, RealFourierNode)
            insizes[subscripts] = 2 * (insizes[subscripts] - 1)
        end
        outranges = [inranges..., chromarange]
        outsizes = [insizes ; layer.bank.spec.n_filters_per_octave]
        for j in layer.bank.behavior.j_range
            ψ_log2_sampling = layer.bank.behavior.ψ_log2_samplings[1+j]
            for subscript in subscripts
                outsizes[subscript] =
                    insizes[subscript] >> (-ψ_log2_sampling)
                inrange = inranges[subscript].second
                outranges[subscript] = outranges[subscript].first =>
                    (inrange.start):(1<<-ψ_log2_sampling):(inrange.stop)
            end
            outdata = zeros(eltype(innode.data), tuple(outsizes...))
            inds = [fill(Colon(), ndims(innode.data)) ; 0]
            for χ in 0:(layer.bank.spec.n_filters_per_octave-1)
                ψ = layer.bank.ψs[1, 1+χ, 1+j]
                inds[end] = 1 + χ
                transform!(view(outdata, inds...), ψ, innode, subscripts...)
            end
            outpath = Path(octavekey => j)
            outnodes[outpath] = Node(outdata, tuple(outranges...))
        end
    end
    WaveletLayerState([ScatteredBlob(outnodes)], layer)
end
