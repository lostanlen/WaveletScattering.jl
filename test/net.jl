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

import WaveletScattering: AbstractNode, get_bandwidth, get_centerfrequency,
    prepend, transform!, scattering_order
