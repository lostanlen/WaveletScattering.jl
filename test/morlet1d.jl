using Base.Test
# morlet1d.jl
import WaveletScattering: Morlet1DSpec, default_motherfreqency,
    fourierwavelet, gauss
# spec.jl
import WaveletScattering: bandwidths, centerfrequencies,
    default_ɛ, gauss, qualityfactors, scales, uncertainty
# meta.jl
import WaveletScattering: NonOrientedMeta

numerictypes = [Float16, Float32, Float64,
                Complex{Float16}, Complex{Float32}, Complex{Float64}]

# Morlet1DSpec default options
for T in numerictypes
  # ordinary defaults, user-specified nOctaves
  spec = Morlet1DSpec(T, nOctaves=8)
  @test_approx_eq spec.ɛ default_ɛ(T)
  @test spec.log2_size == (15,)
  @test_approx_eq spec.max_qualityfactor 1.0
  @test_approx_eq spec.max_scale Inf
  @test spec.nFilters_per_octave == 1
  @test spec.nOctaves == 8
  # nFilters_per_octave defaults to max_qualityfactor when it is provided
  spec = Morlet1DSpec(max_qualityfactor=8)
  @test spec.nFilters_per_octave == 8
  @test spec.nOctaves == 10
  # max_qualityfactor defaults to nFilters_per_octave when it is provided
  spec = Morlet1DSpec(nFilters_per_octave=12)
  @test_approx_eq spec.max_qualityfactor 12.0
  @test spec.nOctaves == 9
end

# Zero-argument constructor
spec = Morlet1DSpec()
@test spec.signaltype == Float32
@test spec.nOctaves == spec.log2_size[1] - 3

# default_motherfrequency
nfos = [1, 2, 4, 8, 12, 16, 24, 32]
for nfo in nfos
    spec = Morlet1DSpec(nFilters_per_octave=nfo)
    ξs = centerfrequencies(spec)
    if nfo==1
        @test_approx_eq ξs[1] 0.39
    else
        @test_approx_eq (ξs[1]-ξs[2]) (1.0 - 2*ξs[1])
    end
end

# fourierwavelet
spec = Morlet1DSpec(nFilters_per_octave = 16)
γs, χs, js = gammas(spec), chromas(spec), octaves(spec)
ξs, qs = centerfrequencies(spec), qualityfactors(spec)
scs, bws = scales(spec), bandwidths(spec)
metas = [NonOrientedMeta(γs[i], χs[i], bws[i], ξs[i], js[i], qs[i], scs[i])
         for i in eachindex(γs)]
for meta in metas
    ψ = fourierwavelet(meta, spec)
    if isa(ψ, Analytic1DFilter)
        @test all(abs(ψ.pos) .> spec.ɛ)
    elseif isa(ψ, VanishingWithMidpoint1DFilter)
        @test all(abs(ψ.an.pos) .> spec.ɛ)
        @test all(abs(ψ.coan.neg) .> spec.ɛ)
        @test abs(ψ.midpoint) > spec.ɛ
    end
end

# gauss
@test_approx_eq gauss(0.0, 1.0) 1.0
for ω in 1.0:10.0
    for den in logspace(0, 3, 4)
        g = gauss(ω, den)
        @test g >= 0.0
        @test_approx_eq g gauss(-ω, den)
        @test_approx_eq sqrt(-log(g) * den) ω
    end
end
