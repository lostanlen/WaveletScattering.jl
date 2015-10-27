# layer.jl
import WaveletScattering: PathKey

# node
subscripts = (PathKey(:time), PathKey(:chunk))
data = rand(Float32, 32768, 128)
node = AbstractFourierNode(data, 1, subscripts)
allocated = @allocated fft!(node)
@test allocated < 4000
