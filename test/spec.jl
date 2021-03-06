using Base.Test
# spec.jl
import WaveletScattering: checkspec_super, default_ɛ,
    default_max_qualityfactor, default_motherfrequency,
    default_n_filters_per_octave, default_n_octaves
# meta.jl
import WaveletScattering: get_γ, get_scale, ΦMeta, ΨMeta1D
# morlet1d.jl
import WaveletScattering: uncertainty
# spec.jl
import WaveletScattering: AbstractSpec
# spec1d.jl
import WaveletScattering: Spec1D
# waveletclass.jl
import WaveletScattering: MexicanHat, Morlet, RedundantWaveletClass

# default_motherfrequency
immutable NullWaveletClass <: RedundantWaveletClass end
class = NullWaveletClass()
for nfo in [1, 2, 4, 8, 12, 24, 32]
    ξ = default_motherfrequency(class, nfo)
    @test 2ξ ≈ (ξ*2^(-1/nfo) + (1-ξ))
end

# default_n_octaves
numerictypes = [Float16, Float32, Float64]
nfos = [1, 2, 4, 8, 12, 24, 32]
for T in numerictypes, nfo in nfos, max_q in nfos[nfos.<=nfo],
    log2_s in (7+ceil(Int, log2(nfo)):18), max_s in [max_q*exp2.(5:14); Inf]
    spec = Spec1D(
        log2_size = log2_s,
        max_qualityfactor = max_q,
        max_scale = max_s,
        n_filters_per_octave = nfo,
        signaltype = T)
    siglength = 1 << log2_s
    if max_s > siglength
        min_centerfrequency = uncertainty(spec) / siglength * max_q
    else
        min_centerfrequency = uncertainty(spec) / max_s * 1.0
    end
    n_octaves = default_n_octaves(nothing, spec.class, log2_s,
                                Float64(max_q), max_s, spec.motherfrequency,
                                nfo)
    ξs = [ meta.centerfrequency for meta in spec.ψmetas ]
    n_octaves_a = floor(Int, log2(ξs[1] / min_centerfrequency))
    n_octaves_b = log2_s - 1 - ceil(Int, log2(nfo))
    @test n_octaves == min(n_octaves_a, n_octaves_b)
end

# checkspec_super
immutable UncheckedSpec <: AbstractSpec
    ɛ::Float64
    log2_size::Int
    ϕmeta::ΦMeta
    ψmetas::Array{ΨMeta1D,3}
    class::RedundantWaveletClass
    max_qualityfactor::Float64
    max_scale::Float64
    motherfrequency::Float64
    n_filters_per_octave::Int
    n_octaves::Int
end
ɛ = 1e-3
log2_size = 13
max_qualityfactor = 12.0
max_scale = 5e3
motherfrequency = 0.45
n_filters_per_octave = 16
n_octaves = 8
ϕmeta = Spec1D().ϕmeta
ψmetas = Spec1D().ψmetas
class = Morlet()
@test_throws ErrorException checkspec_super(UncheckedSpec(-0.1, log2_size,
    ϕmeta, ψmetas, MexicanHat(),
    2.0, max_scale, motherfrequency, n_filters_per_octave,
    n_octaves))
@test_throws ErrorException checkspec_super(UncheckedSpec(-0.1, log2_size,
    ϕmeta, ψmetas, class,
    max_qualityfactor, max_scale, motherfrequency, n_filters_per_octave,
    n_octaves))
@test_throws ErrorException checkspec_super(UncheckedSpec(1.0, log2_size,
    ϕmeta, ψmetas, class,
    max_qualityfactor, max_scale, motherfrequency, n_filters_per_octave,
    n_octaves))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, 1,
    ϕmeta, ψmetas, class,
    max_qualityfactor, max_scale, motherfrequency, n_filters_per_octave,
    n_octaves))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, log2_size,
    ϕmeta, ψmetas, class,
    0.9, max_scale, motherfrequency, n_filters_per_octave, n_octaves))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, log2_size,
    ϕmeta, ψmetas, class,
    max_qualityfactor, max_scale, 0.0, n_filters_per_octave, n_octaves))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, log2_size,
    ϕmeta, ψmetas, class,
    max_qualityfactor, max_scale, 0.51, n_filters_per_octave, n_octaves))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, log2_size,
    ϕmeta, ψmetas, class,
    max_qualityfactor, max_scale, motherfrequency, 0, n_octaves))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, log2_size,
    ϕmeta, ψmetas, class,
    max_qualityfactor, max_scale, motherfrequency, n_filters_per_octave, 0))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, log2_size,
    ϕmeta, ψmetas, class,
    max_qualityfactor, max_scale, motherfrequency, 11, n_octaves))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, log2_size,
    ϕmeta, ψmetas, class,
    16.1, max_scale, motherfrequency, n_filters_per_octave, n_octaves))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, log2_size,
    ϕmeta, ψmetas, class,
    max_qualityfactor, max_scale, motherfrequency, n_filters_per_octave, 14))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, log2_size,
    ϕmeta, ψmetas, class,
    max_qualityfactor, max_scale, motherfrequency, n_filters_per_octave, 12))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, log2_size,
    ϕmeta, ψmetas, class,
    1.0, 1e1, motherfrequency, 1, 12))
@test_throws ErrorException checkspec_super(UncheckedSpec(ɛ, 8,
    ϕmeta, ψmetas, class,
    1.0, Inf, motherfrequency, 1, 1))


# default_ɛ
@test default_ɛ(Float16) ≈ 1e-3
@test default_ɛ(Float32) ≈ 1e-7
@test default_ɛ(Float64) ≈ 1e-15
@test default_ɛ(Complex{Float16}) ≈ 1e-3
@test default_ɛ(Complex{Float32}) ≈ 1e-7
@test default_ɛ(Complex{Float64}) ≈ 1e-15

# default_max_qualityfactor
@test default_max_qualityfactor(8.0, nothing) ≈ 8.0
@test default_max_qualityfactor(nothing, 8) ≈ 8.0
@test default_max_qualityfactor(nothing, nothing) ≈ 1.0

# default_n_filters_per_octave
@test default_n_filters_per_octave(12, nothing) == 12
@test default_n_filters_per_octave(nothing, 7.9) == 8
@test default_n_filters_per_octave(nothing, nothing) == 1

# default_n_octaves (fallback)
type WhateverType end
@test default_n_octaves(5, WhateverType) == 5
