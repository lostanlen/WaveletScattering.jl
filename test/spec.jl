using Base.Test
import WaveletScattering: AbstractSpec, Abstract1DSpec, Abstract2DSpec,
    bandwidths, checkspec, centerfrequencies, chromas, default_ɛ,
    default_max_qualityfactor, default_nFilters_per_octave, default_nOctaves,
    gammas, realtype, octaves, qualityfactors, scales, tune_motherfrequency
import WaveletScattering: Morlet1DSpec

# bandwidths, centerfrequencies, default_nOctaves, qualityfactors, scales
numerictypes = [Float16, Float32, Float64]
nfos = [1, 2, 4, 8, 12, 24, 32]
for T in numerictypes, nfo in nfos, max_q in nfos[nfos.<=nfo],
    log2_s in (4+ceil(Int, log2(nfo)):18), max_s in [max_q*exp2(4:14); Inf]
    machine_precision = max(1e-10, default_ɛ(T))
    spec = Morlet1DSpec(T, nFilters_per_octave=nfo, max_qualityfactor=max_q
                        log2_size=log2_s, max_scale=max_s)
    bws = bandwidths(spec)
    ξs = centerfrequencies(spec)
    qs = qualityfactors(spec)
    scs = scales(spec)
    # bandwidths
    @test_approx_eq bws ξs./qs
    # centerfrequencies
    @test_approx_eq ξs[1] spec.mother_frequency
    difflogξs = diff(log2(ξs))
    @test_approx_eq difflogξs (-ones(difflogξs)/spec.nFilters_per_octave)
    @test all(ξs.>0.0)
    @test_approx_eq ξs bws.*qs
    # default_nOctaves
    siglength = 1 << log2_s
    min_centerfrequency = min(ξs)
    nOctaves = default_nOctaves(nOctaves, Morlet1DSpec, log2_s, max_q, max_s,
        spec.motherfrequency, nfo)
    nOctaves_a = floor(Int, log2(ξs[1] / min(ξs)))
    nOctaves_b = log2_s - 1 - ceil(Int, log2(nfo))
    @test nOctaves = min(nOctaves_a, nOctaves_b)
    # qualityfactors
    qs = qualityfactors(spec)
    @test all(qs.>=0.0)
    @test all(qs.<=max_q)
    @test_approx_eq qs ξs./bws
    # scales
    @test all(scs.>0.0)
    @test all(scs[qs.>1.0] .< (max_s+machine_precision))
    @test all(scs .< (exp2(spec.log2_size[1])+machine_precision))
    # uncertainty
    empirical_uncertainty = bws .* scs
    @test all(abs(diff(empirical_uncertainty)) .< 1e-6)
    @test_approx_eq all(abs(uncertainty(spec)-empirical_uncertainty).<1e-6)
end

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
@test default_nOctaves(5, WhateverType) = 5

# realtype
@test realtype(Float32) == Float32
@test realtype(Float64) == Float64
@test realtype(Complex{Float32}) == Float32
@test realtype(Complex{Float64}) == Float64
@test_throws MethodError realtype(ASCIIString)

# gammas, chromas, octaves
immutable Test1DSpec <: Abstract1DSpec
    nFilters_per_octave::Int
    nOctaves::Int
end
spec = Test1DSpec(1, 1)
@test gammas(spec) == [0]
@test chromas(spec) == [0]
@test octaves(spec) == [0]
spec = Test1DSpec(2, 1)
@test gammas(spec) == [0, 1]
@test chromas(spec) == [0, 1]
@test octaves(spec) == [0, 0]
spec = Test1DSpec(1, 2)
@test gammas(spec) == [0, 1]
@test chromas(spec) == [0, 0]
@test octaves(spec) == [0, 1]
spec = Test1DSpec(12, 8)
nWavelets = spec.nFilters_per_octave * spec.nOctaves
@test length(gammas(spec)) == nWavelets
@test length(chromas(spec)) == nWavelets
@test length(octaves(spec)) == nWavelets

# tune_motherfrequency
nfos = [1, 2, 4, 8, 12, 16, 24]
pitchforks = [392, 415, 422, 430, 435, 440, 442, 444, 466]
for nfo in nfos, pitchfork in pitchforks
    tuning_frequency = pitchfork / 44100.0
    ξ = tune_motherfrequency(spectype, nfo, tuning_frequency)
    spec = Test1DSpec(nFilters_per_octave=nfo)
    γs = gammas(spec)
    ωs = ξ * 2^(-γs / spec.nFilters_per_octave)
    @test any(abs(ωs - ξ) < 1e-4)
    max_ξ = default_motherfrequency(spec)
    @test ξ < max_ξ
    @test ξ * 2^(1/spec.nFilters_per_octave) > max_ξ
end
