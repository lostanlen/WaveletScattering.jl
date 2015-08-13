using Base.Test
import WaveletScattering: AbstractFourier1DFilter, Vanishing1DFilter

N = 8
first = -1
last = 3
y = 1 .<< collect(0:(last-first))
ψ = AbstractFourier1DFilter(y, first, last, N)
@test isa(ψ, Vanishing1DFilter)
@test ψ.an == [1]
@test ψ.coan == [4,8,16]
