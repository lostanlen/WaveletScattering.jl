using Base.Test
# fourierlayerstage.jl
import WaveletScattering: FourierLayerState
# fourierlayer.jl
import WaveletScattering: FourierLayer
# inputlayer.jl
import WaveletScattering: InputLayer
# inputlayerstate.jl
import WaveletScattering: InputLayerState


data = rand(Float32, 256, 2)
backend = Mocha.CPUBackend()
inputlayer = InputLayer(
    tops = [:signal], symbols = [:time, :chunk], data = data)
inputlayerstate = InputLayerState(backend, inputlayer)
