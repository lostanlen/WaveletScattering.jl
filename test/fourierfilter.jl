using Base.Test
# fourierfilter.jl
import WaveletScattering: AbstractFourier1DFilter, Analytic1DFilter,
    Coanalytic1DFilter, FullResolution1DFilter, Symmetric1DFilter,
    Vanishing1DFilter, VanishingWithMidpoint1DFilter,
    getindex, littlewoodpaleyadd!, realtype, renormalize!, scalingfunction,
    spin
# meta.jl
import WaveletScattering: NonOrientedMeta, bandwidths, centerfrequencies,
    chromas, gammas, octaves, qualityfactors, scales
# morlet1d.jl
import WaveletScattering: Morlet1DSpec, fourierwavelet

# multiplication operator with scalar
# Base.(:*){T<:Number}(ψ::Analytic1DFilter{T}, b::Number)
ψ = Analytic1DFilter(Float32[0.1, 0.3], 2)
ψ2 = ψ * 2.0
@test isa(ψ2, Analytic1DFilter{Float32})
@test_approx_eq ψ2.pos Float32[0.2, 0.6]
@test ψ2.posfirst == 2
# Base.(:*){T<:Number}(ψ::Coanalytic1DFilter{T}, b::Number)
ψ = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
ψ2 = ψ * 2.0
@test isa(ψ2, Coanalytic1DFilter{Float32})
@test_approx_eq ψ2.neg Float32[0.2, 0.6, 0.8]
@test ψ2.neglast == -3
# Base.(:*){T<:Number}(ψ::FullResolution1DFilter{T}, b::Number)
ψ = FullResolution1DFilter(Float32[0.01, 0.1, 0.2, 0.3])
ψ2 = ψ * 2.0
@test isa(ψ2, FullResolution1DFilter{Float32})
@test_approx_eq ψ2.coeff Float32[0.02, 0.2, 0.4, 0.6]
# Base.(:*){T<:Number}(ψ::Symmetric1DFilter{T}, b::Number)
ψ = Symmetric1DFilter(Float32[0.8, 0.3, 0.1], one(Float32))
ψ2 = ψ * 2.0
@test isa(ψ2, Symmetric1DFilter{Float32})
@test_approx_eq ψ2.leg Float32[1.6, 0.6, 0.2]
@test_approx_eq ψ2.zero Float32(2.0)
# Base.(:*){T<:Number}(ψ::Vanishing1DFilter{T}, b::Number)
an = Analytic1DFilter(Float32[0.1, 0.3], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
ψ = Vanishing1DFilter(an, coan)
ψ2 = ψ * 2.0
@test isa(ψ2, Vanishing1DFilter{Float32})
@test_approx_eq ψ2.an.pos Float32[0.2, 0.6]
@test ψ2.an.posfirst == 2
@test_approx_eq ψ2.coan.neg Float32[0.2, 0.6, 0.8]
@test ψ2.coan.neglast == -3
# Base.(:*){T<:Number}(ψ::VanishingWithMidpoint1DFilter{T}, b::Number))
an = Analytic1DFilter(Float32[0.1, 0.3], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
midpoint = Float32(0.5)
ψ = VanishingWithMidpoint1DFilter(an, coan, midpoint)
ψ2 = ψ * 2.0
@test isa(ψ2, VanishingWithMidpoint1DFilter{Float32})
@test_approx_eq ψ2.an.pos Float32[0.2, 0.6]
@test ψ2.an.posfirst == 2
@test_approx_eq ψ2.coan.neg Float32[0.2, 0.6, 0.8]
@test ψ2.coan.neglast == -3
@test_approx_eq ψ2.midpoint Float32(1.0)

# getindex
# getindex{T}(ψ::Analytic1DFilter{T}, i::Integer)
ψ = Analytic1DFilter(Float32[0.1, 0.3], 2)
@test Float32[ψ[ω] for ω in 1:4] == Float32[0.0, 0.1, 0.3, 0.0]
# getindex{T}(ψ::Analytic1DFilter{T}, I::UnitRange{Int64})
@test ψ[-1:1] == Float32[0.0 ; 0.0 ; 0.0]
@test ψ[1:4] == Float32[0.0 ; 0.1 ; 0.3; 0.0]
@test ψ[2:3] == Float32[0.1 ; 0.3]
@test ψ[1:3] == Float32[0.0 ; 0.1 ; 0.3]
@test ψ[2:4] == Float32[0.1 ; 0.3 ; 0.0]
@test ψ[5:7] == Float32[0.0 ; 0.0 ; 0.0]
# getindex{T}(ψ::Coanalytic1DFilter{T}, i::Integer)
ψ = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
@test Float32[ψ[ω] for ω in -6:-2] == Float32[0.0, 0.1, 0.3, 0.4, 0.0]
# getindex{T}(ψ::Coanalytic1DFilter{T}, I::UnitRange{Int64})
@test ψ[-3:-1] == Float32[0.4 ; 0.0 ; 0.0]
@test ψ[-4:-3] == Float32[0.3 ; 0.4]
@test ψ[-5:-4] == Float32[0.1 ; 0.3]
@test ψ[-7:-5] == Float32[0.0 ; 0.0 ; 0.1]
@test ψ[-10:-8] == Float32[0.0 ; 0.0 ; 0.0]
@test ψ[2:3] == Float32[0.0 ; 0.0]
@test ψ[-6:-2] == Float32[0.0 ; 0.1 ; 0.3 ; 0.4 ; 0.0]
# getindex{T}(ψ::FullResolution1DFilter{T}, i::Integer)
ψ = FullResolution1DFilter(Float32[0.01, 0.1, 0.2, 0.3])
@test Float32[ψ[ω] for ω in -3:2] == Float32[0.0, 0.2, 0.3, 0.01, 0.1, 0.0]
# getindex{T}(ψ::FullResolution1DFilter{T}, I::UnitRange{Int64})
@test ψ[-8:-6] == Float32[0.0, 0.0, 0.0]
@test ψ[-3:0] == Float32[0.0, 0.2, 0.3, 0.01]
@test ψ[-1:1] == Float32[0.3, 0.01, 0.1]
@test ψ[0:3] == Float32[0.01, 0.1, 0.0, 0.0]
@test ψ[-3:2] == Float32[0.0, 0.2, 0.3, 0.01, 0.1, 0.0]
@test ψ[7:7] == Float32[0.0]
# getindex{T}(ψ::Vanishing1DFilter{T}, i::Integer)
an = Analytic1DFilter(Float32[0.1, 0.3], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
ψ = Vanishing1DFilter(an, coan)
@test [ψ[ω] for ω in -6:4] ==
    Float32[0.0, 0.1, 0.3, 0.4, 0.0, 0.0, 0.0, 0.0, 0.1, 0.3, 0.0]
# getindex{T}(ψ::Vanishing1DFilter{T}, I::UnitRange{Int64})
@test ψ[-6:4] == Float32[0.0, 0.1, 0.3, 0.4, 0.0, 0.0, 0.0, 0.0, 0.1, 0.3, 0.0]
# getindex{T}(ψ::VanishingWithMidpoint1DFilter{T}, I::Integer)
an = Analytic1DFilter(Float32[0.1, 0.3, 0.4], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -2)
midpoint = Float32(0.5)
ψ = VanishingWithMidpoint1DFilter(an, coan, midpoint)
@test Float32[ψ[ω] for ω in -6:5] ==
    Float32[0.0, 0.5, 0.1, 0.3, 0.4, 0.0, 0.0, 0.0, 0.1, 0.3, 0.4, 0.0]
# getindex{T}(ψ::VanishingWithMidpoint1DFilter{T}, I::UnitRange{Int64})
@test ψ[-6:5] ==
    Float32[0.0, 0.5, 0.1, 0.3, 0.4, 0.0, 0.0, 0.0, 0.1, 0.3, 0.4, 0.0]

# littlewoodpaleyadd!
# littlewoodpaleyadd!(lp::Vector, ψ::Analytic1DFilter)
lp = zeros(Float32, 8)
ψ = Analytic1DFilter(Float32[0.1, 0.3], 2)
littlewoodpaleyadd!(lp, ψ); lp = zeros(Float32, 8) # warmup
allocatedmemory = @allocated littlewoodpaleyadd!(lp, ψ)
@test allocatedmemory <= 1e3 # on some machines (e.g. Travis's Linux) it is >0
@test_approx_eq lp Float32[0.0, 0.0, 0.01, 0.09, 0.0, 0.0, 0.0, 0.0]
# littlewoodpaleyadd!(lp::Vector, ψ::Coanalytic1DFilter)
lp = zeros(Float32, 8)
ψ = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
littlewoodpaleyadd!(lp, ψ); lp = zeros(Float32, 8) # warmup
allocatedmemory = @allocated littlewoodpaleyadd!(lp, ψ)
@test allocatedmemory <= 1e3 # on some machines (e.g. Travis's Linux) it is >0
@test_approx_eq lp Float32[0.0, 0.0, 0.0, 0.01, 0.09, 0.16, 0.0, 0.0]
# littlewoodpaleyadd!(lp::Vector, ψ::Symmetric1DFilter)
lp = zeros(Float32, 8)
ψ = Symmetric1DFilter(Float32[0.8, 0.3, 0.1], one(Float32))
littlewoodpaleyadd!(lp, ψ); lp = zeros(Float32, 8) # warmup
allocatedmemory = @allocated littlewoodpaleyadd!(lp, ψ)
@test allocatedmemory <= 1e3 # on some machines (e.g. Travis's Linux) it is >0
@test_approx_eq lp Float32[1.0, 0.64, 0.09, 0.01, 0.0, 0.01, 0.09, 0.64]
# littlewoodpaleyadd!(lp::Vector, ψ::Vanishing1DFilter)
lp = zeros(Float32, 8)
an = Analytic1DFilter(Float32[0.1, 0.3], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
ψ = Vanishing1DFilter(an, coan)
littlewoodpaleyadd!(lp, ψ); lp = zeros(Float32, 8) # warmup
allocatedmemory = @allocated littlewoodpaleyadd!(lp, ψ)
@test allocatedmemory <= 1e3 # on some machines (e.g. Travis's Linux) it is >0
@test_approx_eq lp [0.0, 0.0, 0.01, 0.09, 0.0, 0.01, 0.09, 0.16]
# littlewoodpaleyadd!(lp::Vector, ψ::VanishingWithMidpoint1DFilter)
lp = zeros(Float32, 8)
an = Analytic1DFilter(Float32[0.1, 0.3], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
midpoint = Float32(0.5)
ψ = VanishingWithMidpoint1DFilter(an, coan, midpoint)
littlewoodpaleyadd!(lp, ψ); lp = zeros(Float32, 8) # warmup
allocatedmemory = @allocated littlewoodpaleyadd!(lp, ψ)
@test allocatedmemory <= 1e3 # on some machines (e.g. Travis's Linux) it is >0

# Base.maximum
# Base.maximum(ψ::Analytic1DFilter)
@test_approx_eq maximum(Analytic1DFilter([0.1, -0.3], 2)) 0.3
# Base.maximum(ψ::Coanalytic1DFilter)
@test_approx_eq maximum(Coanalytic1DFilter([0.1, 0.3, 0.4*im], -3)) 0.4
# Base.maximum(ψ::FullResolution1DFilter)
@test_approx_eq maximum(FullResolution1DFilter([0.01, 0.1*im, 0.2, -0.3])) 0.3
# Base.maximum(ψ::Symmetric1DFilter)
@test_approx_eq maximum(Symmetric1DFilter([0.2, 0.1], 1.0)) 1.0
# Base.maximum(ψ::Vanishing1DFilter)
an = Analytic1DFilter(Float32[0.1, 0.3], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
@test_approx_eq maximum(Vanishing1DFilter(an, coan)) Float32(0.4)
# Base.maximum(ψ::VanishingWithMidpoint1DFilter)
an = Analytic1DFilter(Float32[0.1, 0.3], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
midpoint = Float32(0.5)
ψ = VanishingWithMidpoint1DFilter(an, coan, midpoint)
@test_approx_eq maximum(ψ) Float32(0.5)

# renormalize!
# case Q=1
spec = Morlet1DSpec()
γs, χs, js = gammas(spec), chromas(spec), octaves(spec)
ξs, qs = centerfrequencies(spec), qualityfactors(spec)
scs, bws = scales(spec), bandwidths(spec)
@inbounds metas = [
    NonOrientedMeta(γs[i], χs[i], bws[i], ξs[i], js[i], qs[i], scs[i])
    for i in eachindex(γs) ]
@inbounds ψs = AbstractFourier1DFilter{spec.signaltype}[
    fourierwavelet(meta, spec) for meta in metas]
ϕ = scalingfunction(spec)
lp = renormalize!(ψs, ϕ, metas, spec)
@test all(lp.< 1.0001)
N = 1 << spec.log2_size[1]
firstω = round(Int, N * ξs[end])
lastω = round(Int, N * ξs[1])
@test all(lp[1+(firstω:lastω)] .> 0.5)
# case Q>1, max_s = Inf
spec = Morlet1DSpec(nFilters_per_octave=8)
γs, χs, js = gammas(spec), chromas(spec), octaves(spec)
ξs, qs = centerfrequencies(spec), qualityfactors(spec)
scs, bws = scales(spec), bandwidths(spec)
@inbounds metas = [
    NonOrientedMeta(γs[i], χs[i], bws[i], ξs[i], js[i], qs[i], scs[i])
    for i in eachindex(γs) ]
@inbounds ψs = AbstractFourier1DFilter{spec.signaltype}[
    fourierwavelet(meta, spec) for meta in metas]
ϕ = scalingfunction(spec)
lp = renormalize!(ψs, ϕ, metas, spec)
@test all(lp.< 1.001)
N = 1 << spec.log2_size[1]
firstω = round(Int, N * ξs[end])
lastω = round(Int, N * ξs[1])
@test all(lp[1+(firstω:lastω)] .> 0.5)
# case Q>1, max_s < Inf
spec = Morlet1DSpec(nFilters_per_octave=8, max_scale=4410)
γs, χs, js = gammas(spec), chromas(spec), octaves(spec)
ξs, qs = centerfrequencies(spec), qualityfactors(spec)
scs, bws = scales(spec), bandwidths(spec)
@inbounds metas = [
    NonOrientedMeta(γs[i], χs[i], bws[i], ξs[i], js[i], qs[i], scs[i])
    for i in eachindex(γs) ]
@inbounds ψs = AbstractFourier1DFilter{spec.signaltype}[
    fourierwavelet(meta, spec) for meta in metas]
ϕ = scalingfunction(spec)
lp = renormalize!(ψs, ϕ, metas, spec)
@test all(lp.< 1.001)
N = 1 << spec.log2_size[1]
firstω = round(Int, N * ξs[end])
lastω = round(Int, N * ξs[1])
@test all(lp[1+(firstω:lastω)] .> 0.5)

# spin
# spin(::Analytic1DFilter)
ψ = Analytic1DFilter(Float32[0.1, 0.3], 2)
ψspinned = spin(ψ)
@test isa(ψspinned, Coanalytic1DFilter{Float32})
@test_approx_eq ψspinned.neg Float32[0.3, 0.1]
@test ψspinned.neglast == -2
# spin(::Coanalytic1DFilter)
ψ = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
ψspinned = spin(ψ)
@test isa(ψspinned, Analytic1DFilter{Float32})
@test_approx_eq ψspinned.pos Float32[0.4, 0.3, 0.1]
@test ψspinned.posfirst == 3
# spin(::FullResolution1DFilter)
ψ = FullResolution1DFilter(Float32[0.1, 0.2, 0.3, 0.4])
ψspinned = spin(ψ)
@test isa(ψspinned, FullResolution1DFilter{Float32})
@test_approx_eq ψspinned.coeff Float32[0.4, 0.3, 0.2, 0.1]
# spin(::Vanishing1DFilter)
an = Analytic1DFilter(Float32[0.1, 0.3], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
ψ = Vanishing1DFilter(an, coan)
ψspinned = spin(ψ)
@test isa(ψspinned, Vanishing1DFilter{Float32})
@test_approx_eq ψspinned.coan.neg Float32[0.3, 0.1]
@test ψspinned.coan.neglast == -2
@test_approx_eq ψspinned.an.pos Float32[0.4, 0.3, 0.1]
@test ψspinned.an.posfirst == 3
# spin(::VanishingWithMidpoint1DFilter)
an = Analytic1DFilter(Float32[0.1, 0.3], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
midpoint = Float32(0.5)
ψ = VanishingWithMidpoint1DFilter(an, coan, midpoint)
ψspinned = spin(ψ)
@test isa(ψspinned, VanishingWithMidpoint1DFilter{Float32})
@test_approx_eq ψspinned.coan.neg Float32[0.3, 0.1]
@test ψspinned.coan.neglast == -2
@test_approx_eq ψspinned.an.pos Float32[0.4, 0.3, 0.1]
@test ψspinned.an.posfirst == 3
@test_approx_eq ψspinned.midpoint Float32(0.5)
