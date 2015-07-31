import WaveletScattering: Morlet1DSpec, realtype

# Morlet1DSpec constructor
for T in [Float16, Float32, Float64, Complex{Float16}, Complex{Float32}, Complex{Float64}]
    RealT = realtype(T)
    ɛ = 1e-5
    log2_length = 15
    max_qualityfactor = 8.0
    max_scale = 10000.0
    nFilters_per_octave = 12
    nOctaves = 8
    spec = Morlet1DSpec{T}(ɛ, T, log2_length, max_qualityfactor, max_scale,
        nFilters_per_octave, nOctaves)
    @test isa(spec.ɛ, RealT)
    @test spec.signaltype == T
    @test isa(spec.max_qualityfactor, RealT)
    @test isa(spec.max_scale, RealT)
end
