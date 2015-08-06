using Base.Test
import WaveletScattering: AbstractSpec, Abstract1DSpec, Abstract2DSpec,
    bandwidths, checkspec, centerfrequencies, chromas, default_ɛ,
    default_max_qualityfactor, default_nFilters_per_octave, default_nOctaves,
    gammas, realtype, octaves, qualityfactors, scales, tune_motherfrequency
import WaveletScattering: Morlet1DSpec

# bandwidths
numerictypes = [Float16, Float32, Float64]
nfos = [1, 2, 4, 8, 12, 24, 32]
for T in numerictypes, nfo in nfos, max_q in nfos[nfos.<=nfo],
    log2_s in (4+ceil(Int, log2(nfo)):18), max_s in [max_q*exp2(4:14); Inf]
    spec = Morlet1DSpec(T, nFilters_per_octave=nfo, max_qualityfactor=max_q,
                        log2_size=log2_s, max_scale=max_s)

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
spectypes = [Morlet1DSpec]
for nfo in nfos, pitchfork in pitchforks, spectype in spectypes
    tuning_frequency = pitchfork / 44100.0
    ξ = tune_motherfrequency(spectype, nfo, tuning_frequency)
    spec = spectype(nFilters_per_octave=nfo, motherfrequency=ξ)
    γs = gammas(spec)
    ωs = ξ * 2^(-γs / spec.nFilters_per_octave)
    @test any(abs(ωs - ξ) < 1e-4)
    max_ξ = default_motherfrequency(spec)
    @test ξ < max_ξ
    @test ξ * 2^(1/spec.nFilters_per_octave) > max_ξ
end
