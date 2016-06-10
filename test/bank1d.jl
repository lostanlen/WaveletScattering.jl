using Base.Test
# bank1d.jl
import WaveletScattering: Bank1D
# blob.jl
import WaveletScattering: ScatteredBlob
# node.jl
import WaveletScattering: Node
# path.jl
import WaveletScattering: Path
# spec1d.jl
import WaveletScattering: Spec1D

W = Bank1D(Spec1D())
@test ndims(W) == 1

x = zeros(Float32, 1 << W.spec.log2_size)
x[1] = 1.0
Wx = W(x)

@test isa(Wx, ScatteredBlob{Node{Complex{Float32},2}})
