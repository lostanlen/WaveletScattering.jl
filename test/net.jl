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
import WaveletScattering: Log1P, Modulus, PointwiseLayer, PointwiseLayerState
# node.jl
import WaveletScattering: RealFourierNode, InvComplexFourierNode
# waveletlayer.jl
import WaveletScattering: WaveletLayer

J = 8
data = map(Float32, randn(2^J, 2))
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

invfourier = FourierLayer(
    bottoms = [:wavelets],
    tops = [:invfourier],
    pathkeys = [PathKey(:time)]
)

modulus = PointwiseLayer(
    bottoms = [:wavelets],
    tops = [:modulus],
    ρ = Modulus())

log1p = PointwiseLayer(
    bottoms = [:modulus],
    tops = [:log1p],
    ρ = Log1P(Float32(1e-2))
)

layers = Mocha.Layer[signal, fourier, wavelets]

Mocha.init(backend)
net = Mocha.Net("network", backend, layers)

layer = wavelets
inputs = net.output_blobs[:wavelets].nodes

@test isa(net, Mocha.Net{Mocha.CPUBackend})
