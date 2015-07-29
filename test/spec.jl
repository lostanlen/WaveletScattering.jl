immutable TestSpec1D <: AbstractSpec1D
    nFilters_per_octave::Int
    nOctaves::Int
end

# abstract subtype testing
@test TestSpec1D <: AbstractSpec

# bankgammas
spec = TestSpec1D(1, 1)
(γs, χs, js) = bankgammas(spec)
@test γs == [0]
@test χs == [0]
@test js == [0]

spec = TestSpec1D(2, 1)
(γs, χs, js) = bankgammas(spec)
@test γs == [0 1]
@test χs == [0 1]
@test js == [0 0]

spec = TestSpec1D(1, 2)
(γs, χs, js) = bankgammas(spec)
@test γs == [0 1]
@test χs == [0 0]
@test js == [0 1]

spec = TestSpec(12, 8)
(γs, χs, js) = bankgammas(spec)
nWavelets = spec.nFilters_per_octave * spec.nOctaves
@test length(γs) == nWavelets
@test length(χs) == nWavelets
@test length(js) == nWavelets
