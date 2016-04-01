using Base.Test

# using Clp
using DataStructures
using JuMP
using MathProgBase

tests = [
    "domain",
    "group",
    "path",
    "meta",
    "spec",
    "filter",
#    "fourierfilter",
    "weighting",
    "behavior1d",
    "bank1d",
    "morlet1d",
    "pointwise",
    "modulus",
    "symbols",
    "layer",
    "inputlayer",
    "pointwiselayer",
    "fourierlayer",
    "waveletlayer"
]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * ", tfile, " ...")
    include(tfile)
end
