using Base.Test
import WaveletScattering: Morlet1DSpec, bandwidths, centerfrequencies,
    default_ɛ, qualityfactors, scales, uncertainty

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

nfos = [1, 2, 4, 8, 12, 16, 24, 32]
for nfo in nfos
    spec = Morlet1DSpec(nFilters_per_octave=nfo)
    # 1. check that the mother center frequency is at the right place
    ξs = centerfrequencies(spec)
    if nfo==1
        @test_approx_eq ξs[1] 0.39
    else
        @test_approx_eq (ξs[1]-ξs[2]) (1.0 - 2*ξs[1])
    end
    # 2. check that log-frequencies are evenly spaced
    difflogξs = diff(log2(ξs))
    @test_approx_eq difflogξs (-ones(difflogξs)/spec.nFilters_per_octave)
    # 3. check that all center frequencies are strictly positive
    @test all(ξs.>0.0)
end
# qualityfactors, scales, bandwidths, and time-frequency tradeoff.
for T in [Float16, Float32, Float64], nfo in nfos,
  max_q in 1:nfo, max_s in [exp2(7:16); Inf]
  @show T, nfo, max_q, max_s
    try
        spec = Morlet1DSpec(T, max_qualityfactor=max_q, max_scale=max_s,
                            nFilters_per_octave=nfo)
    catch
        continue
    end
    machine_precision = max(1e-10, default_ɛ(T))
    bws = bandwidths(spec)
    ξs = centerfrequencies(spec)
    qs = qualityfactors(spec)
    scs = scales(spec)
    @test all(qs.>=1.0)
    @test all(qs.<=max_q)
    @test all(scs.>0.0)
    @test all(scs[qs.>1.0] .< (max_s+machine_precision))
    @test all(scs .< (exp2(spec.log2_size[1])+machine_precision))
    empirical_uncertainty = bws .* scs
    @test all(abs(diff(empirical_uncertainty)) .< 1e-6)
end
