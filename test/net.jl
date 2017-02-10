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

backend = Mocha.CPUBackend()
signal = InputLayer(
    tops = [:signal],
    symbols = [:time, :chunk],
    data = data)

fourier = FourierLayer(
    bottoms = [:signal],
    tops = [:fourier],
    pathkeys = [PathKey(:time)])

wavelets = WaveletLayer(
    bank = Bank1D(Spec1D(log2_size=J, n_filters_per_octave=4)),
    bottoms = [:fourier],
    tops = [:wavelets]
)

invfourier = InvFourierLayer(
    bottoms = [:wavelets],
    tops = [:invfourier],
    pathkeys = [PathKey(:time)]
)

modulus = PointwiseLayer(
    bottoms = [:invfourier],
    tops = [:modulus],
    ρ = Modulus())

fourier2 = FourierLayer(
    bottoms = [:modulus],
    name = "fourier2",
    tops = [:fourier2],
    pathkeys = [PathKey(:time)]
)

scattering = ScatteringLayer(
    bank = Bank1D(Spec1D(log2_size=J, n_filters_per_octave=2)),
    bottoms = [:fourier2, :wavelets],
    name = "scattering",
    tops = [:scattering]
)

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

import WaveletScattering: get_bandwidth, get_centerfrequency

# setup
layer = scattering
inputs = [fourier2, wavelets]
previous_ψmetas = inputs[2].bank.spec.ψmetas
previous_bandwidths = map(get_bandwidth, previous_ψmetas)
ψmetas = layer.bank.spec.ψmetas
centerfrequencies = map(get_centerfrequency, ψmetas)
n_octaves = size(centerfrequencies, 3)
last_j1s = Array{Int}(n_octaves)

for j2 in 0:(n_octaves-1)
    centerfrequency = centerfrequencies[1, end, 1+j2]
    last_j1s[1+j2] = findlast( sibling_mask_factor *
        previous_bandwidths[1, 1, :] .>= centerfrequency) - 1
    
end
