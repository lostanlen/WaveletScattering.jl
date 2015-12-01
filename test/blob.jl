# blob.jl
import WaveletScattering: ScatteredBlob
# bank.jl
import WaveletScattering: Bank1D
# morlet1d.jl
import WaveletScattering: Spec1D
# path.jl
import WaveletScattering: Path, PathKey
# node.jl
import WaveletScattering: RealFourierNode, InvComplexFourierNode

pathkeys = (PathKey(:time), PathKey(:chunk))
data = rand(Float32, 32768, 256)
node = RealFourierNode(data, [1], pathkeys)
@test_approx_eq maximum(abs(imag(node.data[1,:]))) 0.0
blob = ScatteredBlob(Dict(Path() => node), subscripts)
bank = Bank1D(Spec1D(), PathKey(:time))

invnode = InverseFourierNode(node, [1])

data_out = zeros(Complex{Float32}, 32768, 256)
inverse_plan = plan_ifft!(data_out, 1)
node_out = InverseFourierNode(data_out, inverse_plan, node.ranges)
γ = 0
ψ = bank.ψs[1 + γ]
