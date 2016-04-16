using Base.Test
# inputlayerstate.jl
import WaveletScattering: InputLayerState, kthrange, setup
# inputlayer.jl
import WaveletScattering: InputLayer

backend = Mocha.CPUBackend()
inputlayer = InputLayer(
        tops = [:signal],
        symbols = [:time, :chunk],
        data = rand(Float32, 256, 2))
inputlayerstate_constructor =
    InputLayerState(backend, inputlayer)
@test inputlayerstate_constructor.layer == inputlayer
inputlayerstate_setup = Mocha.setup(backend,
    inputlayer, Mocha.Blob[], Mocha.Blob[])
@test inputlayerstate_setup.layer == inputlayer
