using Base.Test
# fourierfilter.jl
import WaveletScattering: AbstractFourier1DFilter, Analytic1DFilter,
    Coanalytic1DFilter, Vanishing1DFilter, VanishingWithMidpoint1DFilter,
    littlewoodpaleyadd!, realtype, renormalize!, scalingfunction!
# meta.jl
import WaveletScattering: NonOrientedMeta, bandwidths, centerfrequencies,
    chromas, gammas, octaves, qualityfactors, scales
# morlet1d.jl
import WaveletScattering: Morlet1DSpec, fourierwavelet

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
@test_approx_eq lp Float32[0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.09, 0.16]
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

# renormalize!
# case Q=1
spec = Morlet1DSpec()
γs, χs, js = gammas(spec), chromas(spec), octaves(spec)
ξs, qs = centerfrequencies(spec), qualityfactors(spec)
scs, bws = scales(spec), bandwidths(spec)
@inbounds metas = [
    NonOrientedMeta(γs[i], χs[i], bws[i], ξs[i], js[i], qs[i], scs[i])
    for i in eachindex(γs)]
@inbounds ψs = AbstractFourier1DFilter{spec.signaltype}[
    fourierwavelet(meta, spec) for meta in metas]
lp = renormalize!(ψs, metas, spec)
@test all(lp.< 1.001)
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
    for i in eachindex(γs)]
@inbounds ψs = AbstractFourier1DFilter{spec.signaltype}[
    fourierwavelet(meta, spec) for meta in metas]
lp = renormalize!(ψs, metas, spec)
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
    for i in eachindex(γs)]
@inbounds ψs = AbstractFourier1DFilter{spec.signaltype}[
    fourierwavelet(meta, spec) for meta in metas]
lp = renormalize!(ψs, metas, spec)
@test all(lp.< 1.001)
N = 1 << spec.log2_size[1]
firstω = round(Int, N * ξs[end])
lastω = round(Int, N * ξs[1])
# @test all(lp[1+(firstω:lastω)] .> 0.5)

# realtype
@test realtype(Float32) == Float32
@test realtype(Float64) == Float64
@test realtype(Complex{Float32}) == Float32
@test realtype(Complex{Float64}) == Float64
@test_throws MethodError realtype(ASCIIString)
