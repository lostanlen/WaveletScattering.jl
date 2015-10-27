# bank.jl
import WaveletScattering: FourierNonOriented1DBank
# morlet1d.jl
import WaveletScattering: Morlet1DSpec
# layer.jl
import WaveletScattering: PathKey
# node.jl
import WaveletScattering: AbstractFourierNode, InverseFourierNode

# node
subscripts = (PathKey(:time), PathKey(:chunk))
data = rand(Float32, 32768, 256)
node = AbstractFourierNode(data, 1, subscripts)
fft!(node) # warm-up
allocated = @allocated fft!(node)
@test allocated == 0

spec = Morlet1DSpec()
bank = FourierNonOriented1DBank(spec)

data_out = zeros(Complex{Float32}, 32768, 256)
inverse_plan = plan_ifft!(data_out, 1)
node_out = InverseFourierNode(data_out, inverse_plan, node.ranges)
γ = 0
ψ = bank.ψs[1 + γ]
