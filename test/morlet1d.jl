# Morlet1D constructor
for signaltype in [Float16, Float32, Float64,
                   Complex{Float16}, Complex{Float32}, Complex{Float64}]
    RealT = realtype(signaltype)
    ɛ = 1e-5
    log2_length = 15
    max_qualityfactor = 8.0
    max_scale = 10000.0
    nFilters_per_octave = 12
    nOctaves = 8
    m = Morlet1DSpec(ɛ, signaltype, log2_length, max_qualityfactor, max_scale,
        nFilters_per_octave, nOctaves)
    @test isa(m.ɛ, RealT)
    @test m.signaltype == signaltype
    @test isa(max_qualityfactor, RealT)
    @test isa(max_scale, RealT)
end
