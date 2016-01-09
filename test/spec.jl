using Base.Test
# spec.jl
import WaveletScattering: checkspec, default_ɛ, default_max_qualityfactor,
    default_motherfrequency, default_nFilters_per_octave, default_nOctaves,
    tune_motherfrequency
# morlet1d.jl
import WaveletScattering: uncertainty
# spec1d.jl
import WaveletScattering: Spec1D

# default_nOctaves
numerictypes = [Float16, Float32, Float64]
nfos = [1, 2, 4, 8, 12, 24, 32]
for T in numerictypes, nfo in nfos, max_q in nfos[nfos.<=nfo],
    log2_s in (7+ceil(Int, log2(nfo)):18), max_s in [max_q*exp2(5:14); Inf]
    spec = Spec1D(
        signaltype = T,
        log2_size = log2_s,
        max_qualityfactor = max_q,
        max_scale = max_s
        nFilters_per_octave = nfo,)
    siglength = 1 << log2_s
    if max_s > siglength
        min_centerfrequency = uncertainty(Morlet1DSpec) / siglength * max_q
    else
        min_centerfrequency = uncertainty(Morlet1DSpec) / max_s * 1.0
    end
    nOctaves = default_nOctaves(nothing, Morlet1DSpec, tuple(log2_s),
                                Float64(max_q), max_s, spec.motherfrequency,
                                nfo)
    ξs = centerfrequencies(spec)
    nOctaves_a = floor(Int, log2(ξs[1] / min_centerfrequency))
    nOctaves_b = log2_s - 1 - ceil(Int, log2(nfo))
    @test nOctaves == min(nOctaves_a, nOctaves_b)
end

# checkspec
immutable UncheckedSpec <: AbstractSpec
    ɛ::Float64
    log2_size::Tuple{Int}
    max_qualityfactor::Float64
    max_scale::Float64
    motherfrequency::Float64
    nFilters_per_octave::Int
    nOctaves::Int
end
ɛ = 1e-3
log2_size = (13,)
max_qualityfactor = 12.0
max_scale = 5e3
motherfrequency = 0.45
nFilters_per_octave = 16
nOctaves = 8
scales(::UncheckedSpec) = [1e3]
@test_throws ErrorException checkspec(UncheckedSpec(-0.1, log2_size,
    max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave,
    nOctaves))
@test_throws ErrorException checkspec(UncheckedSpec(1.0, log2_size,
    max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave,
    nOctaves))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, (1,),
    max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave,
    nOctaves))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, log2_size,
    0.9, max_scale, motherfrequency, nFilters_per_octave, nOctaves))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, log2_size,
    max_qualityfactor, max_scale, 0.0, nFilters_per_octave, nOctaves))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, log2_size,
    max_qualityfactor, max_scale, 0.51, nFilters_per_octave, nOctaves))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, log2_size,
    max_qualityfactor, max_scale, motherfrequency, 0, nOctaves))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, log2_size,
    max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave, 0))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, log2_size,
    max_qualityfactor, max_scale, motherfrequency, 11, nOctaves))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, log2_size,
    16.1, max_scale, motherfrequency, nFilters_per_octave, nOctaves))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, log2_size,
    max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave, 14))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, log2_size,
    max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave, 12))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, log2_size,
    1.0, 1e1, motherfrequency, 1, 12))
@test_throws ErrorException checkspec(UncheckedSpec(ɛ, (9,),
    max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave, 1))

# default_ɛ
@test_approx_eq default_ɛ(Float16) 1e-3
@test_approx_eq default_ɛ(Float32) 1e-7
@test_approx_eq default_ɛ(Float64) 1e-15
@test_approx_eq default_ɛ(Complex{Float16}) 1e-3
@test_approx_eq default_ɛ(Complex{Float32}) 1e-7
@test_approx_eq default_ɛ(Complex{Float64}) 1e-15

# default_max_qualityfactor
@test_approx_eq default_max_qualityfactor(8.0, nothing) 8.0
@test_approx_eq default_max_qualityfactor(nothing, 8) 8.0
@test_approx_eq default_max_qualityfactor(nothing, nothing) 1.0

# default_nFilters_per_octave
@test default_nFilters_per_octave(12, nothing) == 12
@test default_nFilters_per_octave(nothing, 7.9) == 8
@test default_nFilters_per_octave(nothing, nothing) == 1

# default_nOctaves (fallback)
type WhateverType end
@test default_nOctaves(5, WhateverType) == 5

# tune_motherfrequency
immutable TestSpec <: AbstractSpec
    nFilters_per_octave::Int
    nOctaves::Int
end
nfos = [1, 2, 4, 8, 12, 16, 24]
pitchforks = [392, 415, 422, 430, 435, 440, 442, 444, 466]
for nfo in nfos, pitchfork in pitchforks
    tuningfrequency = pitchfork / 44100.0
    ξ = tune_motherfrequency(tuningfrequency, TestSpec, nfo)
    spec = TestSpec(nfo, 15)
    γs = gammas(spec)
    ωs = ξ * exp2(-γs / nfo)
    @test any(abs(ωs - ξ) .< 1e-4)
    max_ξ = default_motherfrequency(TestSpec, nfo)
    @test ξ < max_ξ
    @test ξ * exp2(inv(nfo)) > max_ξ
end
