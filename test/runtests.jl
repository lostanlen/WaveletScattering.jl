using Base.Test

# using Clp
using DataStructures
using JuMP
using MathProgBase

tests = [
    "domain",
    "spec",
    "meta",
    "filter",
#    "fourierfilter",
    "bank",
    "morlet1d",
    "layer",
    "path"
]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * ", tfile, " ...")
    include(tfile)
end
