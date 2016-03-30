using Base.Test

# using Clp
using DataStructures
using JuMP
using MathProgBase

tests = [
    "domain",
    "path",
    "meta",
    "spec",
    "filter",
#    "fourierfilter",
    "bank",
    "morlet1d",
    "symbols",
    "layer"
]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * ", tfile, " ...")
    include(tfile)
end
