using Base.Test

# meta.jl
import WaveletScattering: get_γ, get_θ, get_χ, get_j,
    get_bandwidth, get_centerfrequency, get_qualityfactor, get_scale

# spec.jl
import WaveletScattering: default_ɛ, uncertainty

# spec1d.jl
import WaveletScattering: Spec1D

# bandwidths, centerfrequencies, qualityfactors, scales, uncertainty
numerictypes = [Float16, Float32, Float64]
nfos = [1, 2, 4, 8, 12, 24, 32]
for T in numerictypes, nfo in nfos, max_q in nfos[nfos.<=nfo],
    log2_s in (7+ceil(Int, log2.(nfo)):18), max_s in [max_q*exp2.(5:14); Inf]
    machine_precision = max(1e-10, default_ɛ(T))
    spec = Spec1D(
            signaltype=T,
            n_filters_per_octave=nfo,
            max_qualityfactor=max_q,
            log2_size=log2_s,
            max_scale=max_s)
    γs = map(get_γ, spec.ψmetas)
    θs = map(get_θ, spec.ψmetas)
    χs = map(get_χ, spec.ψmetas)
    js = map(get_j , spec.ψmetas)
    bws = map(get_bandwidth, spec.ψmetas)
    ξs = map(get_centerfrequency, spec.ψmetas)
    qs = map(get_qualityfactor, spec.ψmetas)
    scs = map(get_scale, spec.ψmetas)
    # gammas
    @test eltype(γs) <: Int16
    # chis
    @test eltype(χs) <: Int8
    # js
    @test eltype(js) <: Int8
    @test map(Int16, χs) + nfo * map(Int16, js) == γs
    # thetas
    @test all(θs .== 0)
    @test eltype(θs) <: Int8
    # bandwidths
    @test_approx_eq bws ξs./qs
    # centerfrequencies
    @test_approx_eq ξs[1] spec.motherfrequency
    difflogξs = diff(log2.(ξs[:]))
    @test difflogξs ≈ (-ones(difflogξs)/spec.n_filters_per_octave)
    @test all(ξs.>0.0)
    @test ξs ≈ bws.*qs
    # qualityfactors
    @test all(qs.>=0.0)
    @test all(qs.<=max_q)
    @test qs ≈ ξs./bws
    # scales
    @test all(scs.>0.0)
    @test all(scs[qs.>1.0] .< (max_s+machine_precision))
    @test all(scs .< (exp2.(spec.log2_size)+machine_precision))
    # uncertainty
    empirical_uncertainty = bws .* scs
    @test all(abs.(diff(empirical_uncertainty[:])) .< machine_precision)
    @test all(abs.(uncertainty(spec)-empirical_uncertainty).< machine_precision)
end
