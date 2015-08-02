import WaveletScattering: Morlet1DSpec, realtype

numerictypes = [Float16, Float32, Float64,
    Complex{Float16}, Complex{Float32}, Complex{Float64}]

# Morlet1DSpec constructor
for T in numerictypes
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

# Morlet1DSpec default options
for T in numerictypes
  opts = @options
  RealT = realtype(T)
  spec = Morlet1DSpec(opts)
  @test spec.ɛ == eps(RealT)
  @test spec.log2_length = 15
  @test_approx_eq spec.max_qualityfactor == 1.0
  @test_approx_eq spec.max_scale == Inf
  @test spec.max_nFilters_per_octave == 1
  @test spec.nOctaves == 8
  opts = @options max_qualityfactor = 8.0
  spec = Morlet1DSpec(opts)
  @test spec.max_nFilters_per_octave == 8
  spec = Morlet1DSpec(opts)
  opts = @options nFilters_per_octave = 12
  @test_approx_eq spec.max_qualityfactor == 12.0
end
