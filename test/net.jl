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
# pooling.jl
import WaveletScattering: Pooling

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
    bank = Bank1D(Spec1D(log2_size=J)),
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

pooling = PoolingLayer(
    bottoms = [:modulus],
    tops = [:pooling],
    pooling = mean
)

layers = Mocha.Layer[signal, fourier, wavelets, invfourier, modulus]

Mocha.init(backend)
net = Mocha.Net("network", backend, layers)

@test isa(net, Mocha.Net{Mocha.CPUBackend})

ws = WaveletScattering


# SUMMARY OF KNOWN BUGS:
