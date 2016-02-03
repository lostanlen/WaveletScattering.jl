using Base.Test
# morlet1d.jl
import WaveletScattering: gauss

# gauss
@test_approx_eq gauss(0.0, 1.0) 1.0
for ω in 1.0:10.0
    for den in logspace(0, 3, 4)
        g = gauss(ω, den)
        @test g >= 0.0
        @test_approx_eq g gauss(-ω, den)
        @test_approx_eq sqrt(-log(g) * den) ω
    end
end
