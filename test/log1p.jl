using Base.Test
import WaveletScattering: Log1P

x = [0.0, 1.0, 2.0, 3.0]
data_in = expm1(x)
data_out = zeros(data_in)
ρ = Log1P(1.0)
Base.map!(ρ, data_out, data_in)
@test_approx_eq x data_out
