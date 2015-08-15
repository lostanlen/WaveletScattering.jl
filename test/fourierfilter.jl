using Base.Test
# fourierfilter.jl
import WaveletScattering: AbstractFourier1DFilter, Analytic1DFilter,
    Coanalytic1DFilter, Vanishing1DFilter, VanishingWithMidpoint1DFilter,
    littlewoodpaleyadd!, realtype

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

# element-wise product .*
# Base.(:.*){T}(ψ::Analytic1DFilter{T}, x::Number)
ψin = Analytic1DFilter(Float32[0.1, 0.3, 0.4], 2)
x = 2.0
ψout = ψin .* x
@test ψout.posfirst == ψin.posfirst
@test isa(ψout.pos, Vector{Float32})
@test_approx_eq ψout.pos Float32[0.2, 0.6, 0.8]
# Base.(:.*){T}(ψ::Analytic1DFilter{T}, x::Vector)
ψin = Analytic1DFilter(Float32[0.1, 0.3, 0.4], 2)
x = collect(0.1:0.1:1.6)
ψout = ψin .* x
@test ψout.posfirst == ψin.posfirst
@test isa(ψout.pos, Vector{Float32})
@test_approx_eq ψout.pos Float32[0.03, 0.12, 0.20]
# Base.(:.*){T}(ψ::Coanalytic1DFilter{T}, x::Number)
ψin = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
x = 2.0
ψout = ψin .* x
@test ψout.neglast == ψin.neglast
@test isa(ψout.neg, Vector{Float32})
@test ψout.neg == Float32[0.2, 0.6, 0.8]
# Base.(:.*){T}(ψ::Coanalytic1DFilter{T}, x::Vector)
ψin = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
x = collect(0.1:0.1:1.6)
ψout = ψin .* x
@test ψout.neglast == ψin.neglast
@test isa(ψout.neg, Vector{Float32})
@test_approx_eq ψout.neg Float32[0.12, 0.39, 0.56]
# Base.(:.*){T}(ψ::Vanishing1DFilter{T}, x::Number)
an = Analytic1DFilter(Float32[0.1, 0.3, 0.4], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.2, 0.3], -3)
ψin = Vanishing1DFilter(an, coan)
x = 2.0
ψout = ψin .* x
@test isa(ψout, Vanishing1DFilter{Float32})
@test_approx_eq ψout.an.pos Float32[0.2, 0.6, 0.8]
@test_approx_eq ψout.coan.neg Float32[0.2, 0.4, 0.6]
# Base.(:.*){T}(ψ::Vanishing1DFilter{T}, x::Vector)
an = Analytic1DFilter(Float32[0.1, 0.3, 0.4], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.2, 0.3], -3)
ψin = Vanishing1DFilter(an, coan)
x = collect(0.1:0.1:1.6)
ψout = ψin .* x
@test isa(ψout, Vanishing1DFilter{Float32})
@test_approx_eq ψout.an.pos Float32[0.03, 0.12, 0.20]
@test_approx_eq ψout.coan.neg Float32[0.12, 0.26, 0.42]
# Base.(:.*){T}(ψ::VanishingWithMidpoint1DFilter{T}, x::Number)
an = Analytic1DFilter(Float32[0.1, 0.3], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.2], -2)
midpoint = Float32(0.5)
ψin = VanishingWithMidpoint1DFilter(an, coan, midpoint)
x = 2.0
ψout = ψin .* x
@test isa(ψout, VanishingWithMidpoint1DFilter{Float32})
@test ψout.an.posfirst == ψin.an.posfirst
@test ψout.coan.neglast == ψin.coan.neglast
@test_approx_eq ψout.an.pos Float32[0.2, 0.6]
@test_approx_eq ψout.coan.neg Float32[0.2, 0.4]
@test_approx_eq ψout.midpoint Float32(1.0)
# Base.(:.*){T}(ψ::VanishingWithMidpoint1DFilter{T}, x::Vector)
an = Analytic1DFilter(Float32[0.1, 0.3], 2)
coan = Coanalytic1DFilter(Float32[0.1, 0.2], -2)
midpoint = Float32(0.5)
ψin = VanishingWithMidpoint1DFilter(an, coan, midpoint)
x = collect(0.1:0.1:0.8)
ψout = ψin .* x
@test isa(ψout, VanishingWithMidpoint1DFilter{Float32})
@test ψout.an.posfirst == ψin.an.posfirst
@test ψout.coan.neglast == ψin.coan.neglast
@test_approx_eq ψout.an.pos Float32[0.03, 0.12]
@test_approx_eq ψout.coan.neg Float32[0.06, 0.14]
@test_approx_eq ψout.midpoint Float32(0.25)

# littlewoodpaleyadd!
# littlewoodpaleyadd!(lp::Vector, ψ::Analytic1DFilter)
lp = zeros(Float32, 8)
ψ = Analytic1DFilter(Float32[0.1, 0.3], 2)
littlewoodpaleyadd!(lp, ψ); lp = zeros(Float32, 8) # warmup
allocatedmemory = @allocated littlewoodpaleyadd!(lp, ψ)
@test allocatedmemory == 0
@test_approx_eq lp Float32[0.0, 0.0, 0.01, 0.09, 0.0, 0.0, 0.0, 0.0]
# littlewoodpaleyadd!(lp::Vector, ψ::Coanalytic1DFilter)
lp = zeros(Float32, 8)
ψ = Coanalytic1DFilter(Float32[0.1, 0.3, 0.4], -3)
littlewoodpaleyadd!(lp, ψ); lp = zeros(Float32, 8) # warmup
allocatedmemory = @allocated littlewoodpaleyadd!(lp, ψ)
@test allocatedmemory == 0
@test_approx_eq lp Float32[0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.09, 0.16]

# realtype
@test realtype(Float32) == Float32
@test realtype(Float64) == Float64
@test realtype(Complex{Float32}) == Float32
@test realtype(Complex{Float64}) == Float64
@test_throws MethodError realtype(ASCIIString)
