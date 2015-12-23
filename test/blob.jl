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

data = rand(Float32, 32768, 256)
backend = Mocha.CPUBackend()
signal_layer =
    InputLayer(name = "signal", data = data, symbols = [:time, :chunk])
signal_state = InputLayerState(backend, signal_layer)

fourier_layer = FourierLayer("fourier", [:data], [:fourier],
    [PathKey(:time)], FFTW.ESTIMATE, Inf)
fourier_state = FourierLayerState(backend, signal_layer, [signal_state])
