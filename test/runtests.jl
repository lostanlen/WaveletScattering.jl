using Base.Test

using Clp
using DataStructures
using JuMP
using MathProgBase

tests = [
    "spec",
    "meta",
    "filter",
    "fourierfilter",
    "bank",
    "morlet1d",
    "path"
]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * ", tfile, " ...")
    include(tfile)
end
