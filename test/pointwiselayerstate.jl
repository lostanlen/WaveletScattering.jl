using Base.Test
# pointwiselayerstate.jl
import WaveletScattering: PointwiseLayerState
# inputlayer.jl
import WaveletScattering: InputLayer
# inputlayerstate.jl
import WaveletScattering: InputLayerState
# path.jl
import WaveletScattering: Path
# pointwiselayer.jl
import WaveletScattering: PointwiseLayer

backend = Mocha.CPUBackend()
data = map(Float32, randn(256, 2))
inputlayer = InputLayer(
        tops = [:signal],
        symbols = [:time, :chunk],
        data = data)
inputlayerstate = InputLayerState(backend, inputlayer)

bottoms = [:signal]
tops = [:identity]
pointwiselayer = PointwiseLayer(bottoms = bottoms, tops = tops)
pointwiselayerstate =
    Mocha.setup(backend, pointwiselayer, inputlayerstate.blobs, Mocha.Blob[])

@test length(pointwiselayerstate.blobs) == 1
pointwisedata = pointwiselayerstate.blobs[1].nodes[Path()].data
@test_approx_eq pointwisedata data
