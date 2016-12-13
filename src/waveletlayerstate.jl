# WaveletLayerState
immutable WaveletLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    blobs::Vector{B}
    blobs_diff::Vector{B}
    layer::WaveletLayer
end

function forward!(backend::Mocha.CPUBackend, state::PointwiseLayerState,
                  ρ::AbstractPointwise, inputs::Vector)
    @inbounds for idblob in eachindex(inputs)
        map!(ρ, state.blobs[idblob], inputs[idblob])
    end
end

"""diffs and inputs must contain `ScatteredBlob`'s."""
function Mocha.setup(
        backend::Mocha.Backend,
        layer::WaveletLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    blobs = Vector{ScatteredBlob}(length(inputs))
    pathkey = layer.bank.behavior.pathkey
    chromakey = prepend(:χ, layer.bank.behavior.pathkey)
    chromarange = chromakey => 0:1:(layer.bank.spec.nFilters_per_octave-1)
    octavekey = prepend(:j, layer.bank.behavior.pathkey)
    for idblob in eachindex(inputs)
        innodes = inputs[idblob].nodes
        outnodes = DataStructures.SortedDict{Path,AbstractFourierNode,
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
            outranges = (collect(inranges)..., chromarange)
            outsizes = [insizes ; layer.bank.spec.nFilters_per_octave]
            for j in layer.bank.behavior.j_range
                ψ_log2_sampling = layer.bank.behavior.ψ_log2_samplings[1+j]
                for subscript in subscripts
                    outsizes[subscript] =
                        insizes[subscript] >> (-ψ_log2_sampling)
                    # TODO update outranges
                end
                outdata = zeros(eltype(innode.data), tuple(outsizes...))
                inds = [fill(Colon(), ndims(innode.data)) ; 0]
                for χ in 0:(layer.bank.spec.nFilters_per_octave-1)
                    ψ = layer.bank.ψs[1, 1+χ, 1+j]
                    inds[end] = 1 + χ
                    transform!(view(outdata, inds...), ψ, innode, subscripts)
                end
                outpath = Path(octavekey => j)
                # should be InvComplexFourierNode instead of Node
                #outnodes[outpath] = Node(outdata, outranges)
            end
        end
        blobs[idblob] = ScatteredBlob(outnodes)
    end
    blobs_diff = blobs
    WaveletLayerState(blobs, blobs_diff, layer)
end
