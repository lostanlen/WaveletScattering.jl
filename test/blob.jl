# layer.jl
import WaveletScattering: InplaceFourierNode, PathKey

subscripts = (PathKey(:time), PathKey(:chunk))
x = rand(Float32, 1024, 16)
node = InplaceFourierNode(x, 1, subscripts)
