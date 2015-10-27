# layer.jl
import WaveletScattering: PathKey
# node.jl
import WaveletScattering.AbstractFourierNode

# node
subscripts = (PathKey(:time), PathKey(:chunk))
data = rand(Float32, 32768, 256)
node = AbstractFourierNode(data, 1, subscripts)
fft!(node) # warm-up
allocated = @allocated fft!(node)
@test allocated == 0
