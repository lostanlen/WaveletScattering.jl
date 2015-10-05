using Base.Test

# meta.jl
import WaveletScattering: bandwidths, centerfrequencies, chromas, gammas,
    octaves, qualityfactors, scales, uncertainty
# spec.jl
import WaveletScattering: default_ɛ, AbstractSpec

# gammas, chromas, octaves
immutable TestSpec <: AbstractSpec
    nFilters_per_octave::Int
    nOctaves::Int
end
spec = TestSpec(1, 1)
@test gammas(spec) == [0]
@test chromas(spec) == [0]
@test octaves(spec) == [0]
spec = TestSpec(2, 1)
@test gammas(spec) == [0, 1]
@test chromas(spec) == [0, 1]
@test octaves(spec) == [0, 0]
spec = TestSpec(1, 2)
@test gammas(spec) == [0, 1]
@test chromas(spec) == [0, 0]
@test octaves(spec) == [0, 1]
spec = TestSpec(12, 8)
nWavelets = spec.nFilters_per_octave * spec.nOctaves
@test length(gammas(spec)) == nWavelets
@test length(chromas(spec)) == nWavelets
@test length(octaves(spec)) == nWavelets

# bandwidths, centerfrequencies, qualityfactors, scales, uncertainty
numerictypes = [Float16, Float32, Float64]
nfos = [1, 2, 4, 8, 12, 24, 32]
for T in numerictypes, nfo in nfos, max_q in nfos[nfos.<=nfo],
    log2_s in (7+ceil(Int, log2(nfo)):18), max_s in [max_q*exp2(5:14); Inf]
    machine_precision = max(1e-10, default_ɛ(T))
    spec = Morlet1DSpec(T, nFilters_per_octave=nfo, max_qualityfactor=max_q,
                        log2_size=log2_s, max_scale=max_s)
    bws = bandwidths(spec)
    ξs = centerfrequencies(spec)
    qs = qualityfactors(spec)
    scs = scales(spec)
    # bandwidths
    @test_approx_eq bws ξs./qs
    # centerfrequencies
    @test_approx_eq ξs[1] spec.motherfrequency
    difflogξs = diff(log2(ξs))
    @test_approx_eq difflogξs (-ones(difflogξs)/spec.nFilters_per_octave)
    @test all(ξs.>0.0)
    @test_approx_eq ξs bws.*qs
    # qualityfactors
    qs = qualityfactors(spec)
    @test all(qs.>=0.0)
    @test all(qs.<=max_q)
    @test_approx_eq qs ξs./bws
    # scales
    @test all(scs.>0.0)
    @test all(scs[qs.>1.0] .< (max_s+machine_precision))
    @test all(scs .< (exp2(spec.log2_size[1])+machine_precision))
    # uncertainty
    empirical_uncertainty = bws .* scs
    @test all(abs(diff(empirical_uncertainty)) .< machine_precision)
    @test all(abs(uncertainty(spec)-empirical_uncertainty).< machine_precision)
end
