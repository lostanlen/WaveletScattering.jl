# blob.jl
import WaveletScattering: ScatteredBlob
# bank.jl
import WaveletScattering: Bank1D
# inputlayer.jl
import WaveletScattering: InputLayer, InputLayerState
# fourierlayer.jl
import WaveletScattering: FourierLayer, FourierLayerState
# morlet1d.jl
import WaveletScattering: Spec1D
# path.jl
import WaveletScattering: Path, PathKey
# node.jl
import WaveletScattering: RealFourierNode, InvComplexFourierNode
# waveletlayer.jl
import WaveletScattering: WaveletLayer
# invfourierlayer.jl
import WaveletScattering: InvFourierLayer
# pointwise.jl
import WaveletScattering: Modulus, PointwiseLayer, PointwiseLayerState
# poolinglayer.jl
import WaveletScattering: PoolingLayer
# scatteringlayer.jl
import WaveletScattering: ScatteringLayer

J = 8
data = zeros(Float32, 2^J, 2)
data[1] = 1.0f0
bank1 = Bank1D(Spec1D(log2_size=J, n_filters_per_octave=1))
bank2 = Bank1D(Spec1D(log2_size=J, n_filters_per_octave=1))

backend = Mocha.CPUBackend()
signal = InputLayer(
    data = data,
    tops = [:signal],
    symbols = [:time, :chunk])

fourier = FourierLayer(
    bottoms = [:signal],
    pathkeys = [PathKey(:time)],
    tops = [:fourier])

wavelets = WaveletLayer(
    bank = bank1,
    bottoms = [:fourier],
    tops = [:wavelets])

invfourier = InvFourierLayer(
    bottoms = [:wavelets],
    pathkeys = [PathKey(:time)],
    tops = [:invfourier])

modulus = PointwiseLayer(
    bottoms = [:invfourier],
    tops = [:modulus],
    ρ = Modulus())

fourier2 = FourierLayer(
    bottoms = [:modulus],
    name = "fourier2",
    pathkeys = [PathKey(:time)],
    tops = [:fourier2])

scattering = ScatteringLayer(
    bank = bank2,
    bottoms = [:fourier2],
    name = "scattering",
    sibling_mask_factor = 1.0,
    previous_bank = bank2,
    tops = [:scattering])

layers = Mocha.Layer[
    signal,
    fourier,
    wavelets,
    invfourier,
    modulus,
    fourier2]

Mocha.init(backend)
net = Mocha.Net("network", backend, layers)

@test isa(net, Mocha.Net{Mocha.CPUBackend})

import WaveletScattering: AbtractNode, get_bandwidth, get_centerfrequency,
    prepend, scattering_order

# setup
layer = scattering
inputs = net.states[6].blobs

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

chromakey = prepend(:χ, layer.bank.behavior.pathkey)
chromarange = chromakey => 0:1:(layer.bank.spec.n_filters_per_octave-1)

outnodes = DataStructures.SortedDict{Path,AbstractNode,
    Base.Order.ForwardOrdering}()

for j2 in layer.bank.behavior.j_range
    centerfrequency = centerfrequencies[1, end, 1+j2]
    bandwidths = layer.sibling_mask_factor * previous_bandwidths[1, 1, :]
    last_j1 = findlast(bandwidths .>= centerfrequency) - 1
    last_j1 = min(last_j1, previous_bank.behavior.j_range.stop)
    log2_sampling = layer.bank.behavior.ψ_log2_samplings[1+j2]
    for inpath in keys(innodes)
        innode = innodes[inpath]
        j1 = inpath.sdict[j1key]
        (j1 < first_j1) && continue
        (j1 > last_j1) && continue
        println(j1, j2)
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
        for χ in 0:(layer.bank.spec.n_filters_per_octave-1)
            ψ = layer.bank.ψs[1, 1+χ, 1+j2]
            println(typeof(ψ))
            inds[end] = 1 + χ
            transform!(view(outdata, inds...), ψ, innode, subscripts...)
        end
        outpath = Path(inpath.sdict..., octavekey => j2)
        outnodes[outpath] = Node(outdata, tuple(outranges...))
    end
end
