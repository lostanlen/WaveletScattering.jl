using Base.Test
# fourierfilter.jl
import WaveletScattering: AbstractFourier1DFilter, Analytic1DFilter,
    Coanalytic1DFilter, Vanishing1DFilter, VanishingWithMidpoint1DFilter,
    littlewoodpaleyadd!, realtype
# meta.jl
import WaveletScattering: NonOrientedMeta
# morlet1d.jl
import WaveletScattering: Morlet1DSpec

# constructors
function test_periodize(y, first, last, log2_length)
    N = 1 << log2_length
    output = zeros(Int, N)
    support = first:last
    for i in eachindex(support)
        n = support[i]
        modn = mod(n, N) + 1
        output[modn] = output[modn] + y[i]
    end
    output[1] = 0
    return output
end

log2_length = 4
N = 1 << log2_length

first = -2
last = 3
y = 1 .<< collect(0:(last-first))
ψ = AbstractFourier1DFilter(y, first, last, log2_length)
@test isa(ψ, Vanishing1DFilter)
periodized_y = zeros(Int, N)
periodized_y[ψ.an.posfirst + (1:length(ψ.an.pos))] = ψ.an.pos
periodized_y[(N+ψ.coan.neglast+3) + ((-length(ψ.coan.neg)):(-1))] = ψ.coan.neg
@test periodized_y == test_periodize(y, first, last, log2_length)

first = -7
last = 15
y = 1 .<< collect(0:(last-first))
ψ = AbstractFourier1DFilter(y, first, last, log2_length)
@test isa(ψ, VanishingWithMidpoint1DFilter)
periodized_y = zeros(Int, N)
periodized_y[ψ.an.posfirst + (1:length(ψ.an.pos))] = ψ.an.pos
periodized_y[(N+ψ.coan.neglast+2) + ((-length(ψ.coan.neg)):(-1))] = ψ.coan.neg
periodized_y[1 + end >> 1] = ψ.midpoint
@test periodized_y == test_periodize(y, first, last, log2_length)

first = -2
last = 9
y = 1 .<< collect(0:(last-first))
ψ = AbstractFourier1DFilter(y, first, last, log2_length)
@test isa(ψ, VanishingWithMidpoint1DFilter)
periodized_y = zeros(Int, N)
periodized_y[ψ.an.posfirst + (1:length(ψ.an.pos))] = ψ.an.pos
periodized_y[(N+ψ.coan.neglast+2) + ((-length(ψ.coan.neg)):(-1))] = ψ.coan.neg
periodized_y[1 + end >> 1] = ψ.midpoint
@test periodized_y == test_periodize(y, first, last, log2_length)

first = -2
last = 17
y = 1 .<< collect(0:(last-first))
ψ = AbstractFourier1DFilter(y, first, last, log2_length)
@test isa(ψ, VanishingWithMidpoint1DFilter)
periodized_y = zeros(Int, N)
periodized_y[ψ.an.posfirst + (1:length(ψ.an.pos))] = ψ.an.pos
periodized_y[(N+ψ.coan.neglast+2) + ((-length(ψ.coan.neg)):(-1))] = ψ.coan.neg
periodized_y[1 + end >> 1] = ψ.midpoint
@test periodized_y == test_periodize(y, first, last, log2_length)

first = -2
last = 24
y = 1 .<< collect(0:(last-first))
ψ = AbstractFourier1DFilter(y, first, last, log2_length)
@test isa(ψ, VanishingWithMidpoint1DFilter)
periodized_y = zeros(Int, N)
periodized_y[ψ.an.posfirst + (1:length(ψ.an.pos))] = ψ.an.pos
periodized_y[(N+ψ.coan.neglast+2) + ((-length(ψ.coan.neg)):(-1))] = ψ.coan.neg
periodized_y[1 + end >> 1] = ψ.midpoint
@test periodized_y == test_periodize(y, first, last, log2_length)

first = 2
last = 3
y = 1 .<< collect(0:(last-first))
ψ = AbstractFourier1DFilter(y, first, last, log2_length)
@test isa(ψ, Analytic1DFilter)
periodized_y = zeros(Int, N)
periodized_y[ψ.posfirst + (1:length(ψ.pos))] = ψ.pos
@test periodized_y == test_periodize(y, first, last, log2_length)

first = 2
last = 15
y = 1 .<< collect(0:(last-first))
ψ = AbstractFourier1DFilter(y, first, last, log2_length)
@test isa(ψ, VanishingWithMidpoint1DFilter)
periodized_y = zeros(Int, N)
periodized_y[ψ.an.posfirst + (1:length(ψ.an.pos))] = ψ.an.pos
periodized_y[(N+ψ.coan.neglast+2) + ((-length(ψ.coan.neg)):(-1))] = ψ.coan.neg
periodized_y[1 + end >> 1] = ψ.midpoint
@test periodized_y == test_periodize(y, first, last, log2_length)

first = 2
last = 17
y = 1 .<< collect(0:(last-first))
ψ = AbstractFourier1DFilter(y, first, last, log2_length)
@test isa(ψ, VanishingWithMidpoint1DFilter)
periodized_y = zeros(Int, N)
periodized_y[ψ.an.posfirst + (1:length(ψ.an.pos))] = ψ.an.pos
periodized_y[(N+ψ.coan.neglast+2) + ((-length(ψ.coan.neg)):(-1))] = ψ.coan.neg
periodized_y[1 + end >> 1] = ψ.midpoint
@test periodized_y == test_periodize(y, first, last, log2_length)

first = 7
last = 17
y = 1 .<< collect(0:(last-first))
ψ = AbstractFourier1DFilter(y, first, last, log2_length)
@test isa(ψ, VanishingWithMidpoint1DFilter)
periodized_y = zeros(Int, N)
periodized_y[ψ.an.posfirst + (1:length(ψ.an.pos))] = ψ.an.pos
periodized_y[(N+ψ.coan.neglast+2) + ((-length(ψ.coan.neg)):(-1))] = ψ.coan.neg
periodized_y[1 + end >> 1] = ψ.midpoint
@test periodized_y == test_periodize(y, first, last, log2_length)

first = 2
last = 24
y = 1 .<< collect(0:(last-first))
ψ = AbstractFourier1DFilter(y, first, last, log2_length)
@test isa(ψ, VanishingWithMidpoint1DFilter)
periodized_y = zeros(Int, N)
periodized_y[ψ.an.posfirst + (1:length(ψ.an.pos))] = ψ.an.pos
periodized_y[(N+ψ.coan.neglast+2) + ((-length(ψ.coan.neg)):(-1))] = ψ.coan.neg
periodized_y[1 + end >> 1] = ψ.midpoint
@test periodized_y == test_periodize(y, first, last, log2_length)

# right division /
x = 2.0
# Base.(:/){T}(ψ::Analytic1DFilter{T}, x::Number)
an = Analytic1DFilter(Float32[1.2, 0.8, 1.4], 2)
ψout = an / x
@test ψout.posfirst == an.posfirst
@test isa(ψout.pos, Vector{Float32})
@test_approx_eq ψout.pos Float32[0.6, 0.4, 0.7]
# Base.(:/){T}(ψ::Coanalytic1DFilter{T}, x::Number)
coan = Coanalytic1DFilter(Float32[1.2, 0.8, 1.4], -2)
ψout = coan / x
@test ψout.neglast == coan.neglast
@test isa(ψout.neg, Vector{Float32})
@test_approx_eq ψout.neg Float32[0.6, 0.4, 0.7]
# Base.(:/){T}(ψ::Vanishing1DFilter{T}, x::Number)
vanishing = Vanishing1DFilter(an, coan)
ψout = vanishing / x
@test isa(ψout, Vanishing1DFilter)
# Base.(:/){T}(ψ::Vanishing1DFilter{T}, x::Number)
vanishingwithmidpoint = VanishingWithMidpoint1DFilter(an, coan, Float32(0.5))
ψout = vanishingwithmidpoint / x
@test isa(ψout, VanishingWithMidpoint1DFilter)

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
@inbounds ψs = [fourierwavelet(meta, spec) for meta in metas]
lp = renormalize!(ψs, metas, spec)
@test all(lp.<1.0)
N = 1 << log2_length[1]
firstω = round(Int, N * ξs[end])
lastω = round(Int, N * ξs[1])
@test all(lp[1+(firstω:lastω)] .> 0.5)
# case Q>1, max_s = Inf

# case Q>1, max_s < Inf

# realtype
@test realtype(Float32) == Float32
@test realtype(Float64) == Float64
@test realtype(Complex{Float32}) == Float32
@test realtype(Complex{Float64}) == Float64
@test_throws MethodError realtype(ASCIIString)
