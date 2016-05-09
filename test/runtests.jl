using Base.Test

# using Clp
using DataStructures
using JuMP
using MathProgBase

tests = [
    "domain",
    "group",
    "waveletclass",
    "path",
    "meta",
    "spec",
    "spec1d",
    "spec2d",
    "filter",
    "fourierfilter",
    "weighting",
    "behavior1d",
    "bank1d",
    "morlet1d",
    "pointwise",
    "modulus",
    "log1p",
    "symbols",
    "layer",
    "inputlayer",
    "pointwiselayer",
    "fourierlayer",
    "waveletlayer",
    "inputlayerstate",
    "pointwiselayerstate",
    "fourierlayerstate",
    "net"
]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * ", tfile, " ...")
    include(tfile)
end
