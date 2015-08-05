import WaveletScattering: AbstractSpec, Abstract1DSpec, Abstract2DSpec,
    checkspec, chromas, default_ɛ, default_max_qualityfactor,
    default_nFilters_per_octave, default_nOctaves, gammas, realtype, octaves,
    tune

# checkspec
ɛ = eps(Float32)
log2_length = 15
max_qualityfactor = 12.0
max_scale = 1e5
nFilters_per_octave = 12
nOctaves = 8
# normal behavior
@test checkspec(ɛ, log2_length, max_qualityfactor,
                max_scale, nFilters_per_octave, nOctaves)
# error if ɛ is strictly negative
@test_throws ErrorException checkspec(-1.0, log2_length, max_qualityfactor,
                                      max_scale, nFilters_per_octave, nOctaves)
# error if ɛ equal or greater than one
@test_throws ErrorException checkspec(1.0, log2_length, max_qualityfactor,
                                      max_scale, nFilters_per_octave, nOctaves)
# error if log2_length is smaller than 2
@test_throws ErrorException checkspec(ɛ, 1, max_qualityfactor,
                                      max_scale, nFilters_per_octave, nOctaves)
# error if nOctaves is smaller than 1
@test_throws ErrorException checkspec(ɛ, log2_length, max_qualityfactor,
                                      max_scale, nFilters_per_octave, 0)
# error if nOctave >= log2_length
@test_throws ErrorException checkspec(ɛ, 9, max_qualityfactor,
                                      max_scale, 1, 9)
# error if log2_length-nOctaves <= log2(nFilters_per_octave)
@test_throws ErrorException checkspec(ɛ, 12, max_qualityfactor,
                                      max_scale, 32, 9)
# error if max_qualityfactor > nFilters_per_octave
@test_throws ErrorException checkspec(ɛ, log2_length, 8.0,
                                      max_scale, 4, nOctaves)
# error if max_qualityfactor is strictly smaller than 1.0
@test_throws ErrorException checkspec(ɛ, log2_length, 0.9,
                                      max_scale, nFilters_per_octave, nOctaves)
# error if max_scale < (5.0*max_qualityfactor)
@test_throws ErrorException checkspec(ɛ, log2_length, 32.0,
                                      100.0, nFilters_per_octave, nOctaves)
# error if nFilters_per_octave < 1
@test_throws ErrorException checkspec(ɛ, log2_length, max_qualityfactor,
                                      max_scale, 0, nOctaves)

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
