# ScatteringLayerState
immutable ScatteringLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    blobs::Vector{B}
    layer::ScatteringLayer
end

function ScatteringLayerState(
        backend::Mocha.Backend,
        layer::ScatteringLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    innodes = inputs[1].nodes
    firstpath = collect(keys(innodes))[1]
    pathkey = layer.bank.behavior.pathkey
    order = scattering_order(firstpath, pathkey)
    assert(order > 0)
    j1key = prepend((:j, order), pathkey)
    previous_bank = layer.previous_bank
    previous_ψmetas = previous_bank.spec.ψmetas
    previous_bandwidths = map(get_bandwidth, previous_ψmetas)
    previous_log2_samplings = previous_bank.behavior.ψ_log2_samplings
    ψmetas = layer.bank.spec.ψmetas
    centerfrequencies = map(get_centerfrequency, ψmetas)
    n_octaves = size(centerfrequencies, 3)
    first_j1 = layer.previous_bank.behavior.j_range.start
    paths = collect(keys(innodes))
    chromakey = prepend((:χ, 1+order), pathkey)
    n_filters_per_octave = layer.bank.spec.n_filters_per_octave
    chromarange = chromakey => 0:1:(n_filters_per_octave-1)
    octavekey = prepend((:j, 1+order), pathkey)
    bandwidths = layer.sibling_mask_factor * previous_bandwidths[1, 1, :]
    outnodes = DataStructures.SortedDict{Path,AbstractNode,
        Base.Order.ForwardOrdering}()
    for j2 in layer.bank.behavior.j_range
        centerfrequency = centerfrequencies[1, end, 1+j2]
        last_j1 = findlast(bandwidths .>= centerfrequency) - 1
        last_j1 = min(last_j1, previous_bank.behavior.j_range.stop)
        log2_sampling = layer.bank.behavior.ψ_log2_samplings[1+j2]
        for inpath in keys(innodes)
            innode = innodes[inpath]
            j1 = inpath.sdict[j1key]
            (j1 < first_j1) && continue
            (j1 > last_j1) && continue
            inranges = innode.ranges
            inpathkeys = [ pair.first for pair in collect(inranges) ]
            subscripts = find(pathkey .== inpathkeys)
            insizes = collect(size(innode.data))
            if isa(innode, RealFourierNode)
                insizes[subscripts] = 2 * (insizes[subscripts] - 1)
            end
            outranges = [inranges..., chromarange]
            outsizes = [insizes ; layer.bank.spec.n_filters_per_octave]
            previous_log2_sampling = previous_log2_samplings[1+j1]
            log2_resampling = log2_sampling - previous_log2_sampling
            for subscript in subscripts
                outsizes[subscript] = insizes[subscript] >> (-log2_resampling)
                inrange = inranges[subscript].second
                outranges[subscript] = outranges[subscript].first =>
                    (inrange.start):(1<<-log2_resampling):(inrange.stop)
            end
            outdata = zeros(eltype(innode.data), tuple(outsizes...))
            inds = [fill(Colon(), ndims(innode.data)) ; 0]
            for χ2 in 0:(layer.bank.spec.n_filters_per_octave-1)
                ψ = layer.bank.ψs[1, 1+χ2, 1+j2]
                inds[end] = 1 + χ2
                transform!(view(outdata, inds...), ψ, innode, subscripts...)
            end
            outpath = Path(inpath.sdict..., octavekey => j2)
            outnodes[outpath] = Node(outdata, tuple(outranges...))
        end
    end
    ScatteringLayerState([ScatteredBlob(outnodes)], layer)
end

function Mocha.setup(
        backend::Mocha.Backend,
        layer::ScatteringLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    return ScatteringLayerState(backend, layer, inputs, diffs)
end
