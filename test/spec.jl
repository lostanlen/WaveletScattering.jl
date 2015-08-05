import WaveletScattering: AbstractSpec, Abstract1DSpec, Abstract2DSpec,
    checkspec, chromas, default_ɛ, default_max_qualityfactor,
    default_nFilters_per_octave, default_nOctaves, gammas, realtype, octaves,
    tune
import WaveletScattering: Morlet1DSpec

# checkspec
ɛ = eps(Float32)
log2_size = 15
max_qualityfactor = 12.0
max_scale = 1e5
motherfrequency = 0.48
nFilters_per_octave = 12
nOctaves = 8
# normal behavior
@test
    checkspec(ɛ, log2_size, max_qualityfactor,
              max_scale, motherfrequency, nFilters_per_octave, nOctaves)
# error if ɛ < 0.0
@test_throws ErrorException
    checkspec(-0.1, log2_size, max_qualityfactor,
              max_scale, motherfrequency, nFilters_per_octave, nOctaves)
# error if ɛ >= 1.0
@test_throws ErrorException
    checkspec(1.0, log2_size, max_qualityfactor,
              max_scale, motherfrequency, nFilters_per_octave, nOctaves)
# error if log2_size < 2
@test_throws ErrorException
    checkspec(ɛ, 1, max_qualityfactor,
              max_scale, motherfrequency, nFilters_per_octave, nOctaves)
# error if max_qualityfactor < 1.0
@test_throws ErrorException
    checkspec(ɛ, log2_size, 0.9,
              max_scale, motherfrequency, nFilters_per_octave, nOctaves)
# error if motherfrequency <= 0
@test_throws ErrorException
    checkspec(ɛ, log2_size, max_qualityfactor,
              max_scale, 0.0, nFilters_per_octave, nOctaves)
# error if motherfrequency > 0.5
@test_throws ErrorException
    checkspec(ɛ, log2_size, max_qualityfactor,
              max_scale, 0.51, nFilters_per_octave, nOctaves)
# error if nFilters_per_octave < 1
@test_throws ErrorException
    checkspec(ɛ, log2_size, max_qualityfactor,
              max_scale, motherfrequency, 0, nOctaves)
# error if nOctaves < 1
@test_throws ErrorException
    checkspec(ɛ, log2_size, max_qualityfactor,
              max_scale, motherfrequency, nFilters_per_octave, 0)
# error if nOctave >= log2_size
@test_throws ErrorException
    checkspec(ɛ, 9, max_qualityfactor,
              max_scale, motherfrequency, 1, 9)
# error if log2_size-nOctaves <= log2(nFilters_per_octave)
@test_throws ErrorException
    checkspec(ɛ, 12, max_qualityfactor,
              max_scale, motherfrequency, 32, 9)
# error if max_qualityfactor > nFilters_per_octave
@test_throws ErrorException
    checkspec(ɛ, log2_size, 8.0,
              max_scale, motherfrequency, 4, nOctaves)
# error if max_scale < (5.0*max_qualityfactor)
@test_throws ErrorException
    checkspec(ɛ, log2_size, 32.0,
              100.0, motherfrequency, nFilters_per_octave, nOctaves)

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
@test default_nFilters_per_octave(nothing, 12) == 12.0
@test default_nFilters_per_octave(7.9, nothing) == 8
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

# tune
nfos = [1, 2, 4, 8, 12, 16, 24]
pitchforks = [392, 415, 422, 430, 435, 440, 442, 444, 466]
spectypes = [Morlet1DSpec]
for nfo in nfos, pitchfork in pitchforks, spectype in spectypes
    tuning_frequency = pitchfork / 44100.0
    ξ = tune(spectype, nfo, tuning_frequency)
    spec = spectype(nFilters_per_octave=nfo, motherfrequency=ξ)
    γs = gammas(spec)
    ωs = ξ * 2^(-γs / spec.nFilters_per_octave)
    @test any(abs(ωs - ξ) < 1e-4)
    max_ξ = default_motherfrequency(spec)
    @test ξ < max_ξ
    @test ξ * 2^(1/spec.nFilters_per_octave) > max_ξ
end
