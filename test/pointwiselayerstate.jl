using Base.Test
# pointwiselayerstate.jl
import WaveletScattering: PointwiseLayerState
# inputlayer.jl
import WaveletScattering: InputLayer
# inputlayerstate.jl
import WaveletScattering: InputLayerState
# pointwiselayer.jl
import WaveletScattering: PointwiseLayer

backend = Mocha.CPUBackend()
inputlayer = InputLayer(
        tops = [:signal],
        symbols = [:time, :chunk],
        data = map(Float32, randn(256, 2)))
inputlayerstate = InputLayerState(backend, inputlayer)

bottoms = [:signal]
tops = [:identity]
pointwiselayer = PointwiseLayer(bottoms=bottoms, tops=tops)



diffs = Mocha.Blob[]
pointwiselayerstate =
    PointwiseLayerState(backend, layer, inputs, diffs)
