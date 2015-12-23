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

data = rand(Float32, 32768, 256)
backend = Mocha.CPUBackend()
signal_layer = InputLayer(
        top = :signal,
        symbols = [:time, :chunk],
        data = data)
signal_state = InputLayerState(backend, signal_layer)

fourier_layer = FourierLayer(
        bottoms = [:signal],
        tops = [:fourier],
        pathkeys = [PathKey(:time)])
fourier_state = FourierLayerState(backend, fourier_layer, signal_state.blobs)

abs_layer = PointwiseLayer(
        bottoms = [:fourier],
        tops = [:modulus],
        ρ = Modulus())
abs_state = PointwiseLayerState(backend, abs_layer, fourier_state.blobs)
