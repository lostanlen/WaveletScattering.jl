import WaveletScattering: AbstractSpec, Abstract1DSpec, Abstract2DSpec,
                          checkspec, default_epsilon, realtype, specgammas

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
# error if ɛ is strictly negative zero
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
# error if the difference (log2_length-nOctaves) is smaller than 2
@test_throws ErrorException checkspec(ɛ, 9, max_qualityfactor,
    max_scale, 1, 8)
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

# default_epsilon
@test_approx_eq default_epsilon(Float16) 1e-3
@test_approx_eq default_epsilon(Float32) 1e-7
@test_approx_eq default_epsilon(Float64) 1e-15
@test_approx_eq default_epsilon(Complex{Float16}) 1e-3
@test_approx_eq default_epsilon(Complex{Float32}) 1e-7
@test_approx_eq default_epsilon(Complex{Float64}) 1e-15

# realtype
@test realtype(Float32) == Float32
@test realtype(Float64) == Float64
@test realtype(Complex{Float32}) == Float32
@test realtype(Complex{Float64}) == Float64
@test_throws MethodError realtype(ASCIIString)

# specgammas
immutable Test1DSpec <: Abstract1DSpec
    nFilters_per_octave::Int
    nOctaves::Int
end
spec = Test1DSpec(1, 1)
(γs, χs, js) = specgammas(spec)
@test γs == [0]
@test χs == [0]
@test js == [0]

spec = Test1DSpec(2, 1)
(γs, χs, js) = specgammas(spec)
@test γs == [0, 1]
@test χs == [0, 1]
@test js == [0, 0]

spec = Test1DSpec(1, 2)
(γs, χs, js) = specgammas(spec)
@test γs == [0, 1]
@test χs == [0, 0]
@test js == [0, 1]

spec = Test1DSpec(12, 8)
(γs, χs, js) = specgammas(spec)
nWavelets = spec.nFilters_per_octave * spec.nOctaves
@test length(γs) == nWavelets
@test length(χs) == nWavelets
@test length(js) == nWavelets
