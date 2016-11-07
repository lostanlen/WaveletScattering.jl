using Base.Test
# spec.jl
import WaveletScattering: tune_motherfrequency, default_motherfrequency
# meta.jl
import WaveletScattering: get_γ
# spec1d.jl
import WaveletScattering: Spec1D

# tune_motherfrequency
nfos = [1, 2, 4, 8, 12, 16, 24]
pitchforks = [392, 415, 422, 430, 435, 440, 442, 444, 466]
for nfo in nfos, pitchfork in pitchforks
    tuningfrequency = pitchfork / 44100.0
    spec = Spec1D()
    ξ = tune_motherfrequency(tuningfrequency, spec.class, nfo)
    γs = map(get_γ, spec.ψmetas)
    ωs = ξ * exp2.(-γs / nfo)
    @test any(abs.(ωs - ξ) .< 1e-4)
    max_ξ = default_motherfrequency(spec.class, nfo)
    @test ξ < max_ξ
    @test ξ * exp2(inv(nfo)) > max_ξ
end
