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
    "bank",
    "morlet1d",
    "modulus",
    "symbols",
    "layer",
    "fourierlayer"
]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * ", tfile, " ...")
    include(tfile)
end
