import WaveletScattering: AbstractSpec, Abstract1DSpec, Abstract2DSpec,
                          checkspec, realtype, specgammas

immutable Test1DSpec <: Abstract1DSpec
    nFilters_per_octave::Int
    nOctaves::Int
end
immutable Test2DSpec <: Abstract2DSpec
    nFilters_per_octave::Int
    nOctaves::Int
    nOrientations::Int
end

# abstract subtype testing
@test Test1DSpec <: AbstractSpec
@test Test2DSpec <: AbstractSpec

# checkspec
max_qualityfactor = 2.0
nFilters_per_octave = 4
@test checkspec(max_qualityfactor, nFilters_per_octave) == nothing
max_qualityfactor = 0.5
@test_throws ErrorException checkspec(max_qualityfactor, nFilters_per_octave)
max_qualityfactor = 8.0
@test_throws ErrorException checkspec(max_qualityfactor, nFilters_per_octave)

# realtype
@test realtype(Float32) == Float32
@test realtype(Float64) == Float64
@test realtype(Complex{Float32}) == Float32
@test realtype(Complex{Float64}) == Float64
@test_throws MethodError realtype(ASCIIString)

# specgammas
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
