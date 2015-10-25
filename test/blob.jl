# layer.jl
import WaveletScattering: FourierNode, PathKey

# node
subscripts = (PathKey(:time), PathKey(:chunk))
x = rand(Float32, 1024, 16)
node = InplaceFourierNode(x, 1, subscripts)
allocated = @allocated fft!(node)
@test allocated < 4000
