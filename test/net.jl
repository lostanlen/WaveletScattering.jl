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
# pointwise.jl
import WaveletScattering: Modulus, PointwiseLayer, PointwiseLayerState
# node.jl
import WaveletScattering: RealFourierNode, InvComplexFourierNode

data = map(Float32, randn(256, 2))
backend = Mocha.CPUBackend()
signal = InputLayer(
    tops = [:signal],
    symbols = [:time, :chunk],
    data = data)

fourier = FourierLayer(
    bottoms = [:signal],
    tops = [:fourier],
    pathkeys = [PathKey(:time)])

modulus = PointwiseLayer(
    bottoms = [:signal],
    tops = [:modulus],
    ρ = Modulus())

layers = Mocha.Layer[signal]

Mocha.init(backend)
net = Mocha.Net("fourier-modulus", backend, layers)
