using Base.Test
# pointwise.jl
import WaveletScattering: Identity, Log1P
# blob.jl
import WaveletScattering: ScatteredBlob
# inputlayer.jl
import WaveletScattering: InputLayer
# inputlayerstate.jl
import WaveletScattering: InputLayerState, kthrange
# node.jl
import WaveletScattering: Node
# path.jl
import WaveletScattering: Path

# Identity
ρ = Identity()
@test ρ([1.0, 2.0, 3.0]) ≈ [1.0, 2.0, 3.0]

# Log1P
threshold = 0.1
ρ = Log1P(threshold)
@test ρ([0.0]) ≈ [0.0]
@test expm1.(ρ([1.0, 2.0, 3.0])) ≈ [1.0, 2.0, 3.0]*threshold
@test_throws DomainError ρ([- 1.001 / threshold])

# Base.map
layer = InputLayer(
        tops = [:signal],
        symbols = [:time, :chunk],
        data = map(Float32, randn(256, 2)))
ranges = ntuple(k -> kthrange(layer, k), ndims(layer.data))
blob = ScatteredBlob(
    DataStructures.SortedDict((Path() => Node(layer.data, ranges),)))
ρ = Identity()
outblob = ρ(blob)
