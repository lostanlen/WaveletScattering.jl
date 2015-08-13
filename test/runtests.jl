using Base.Test

using DataStructures

tests = [
    "spec",
    "meta",
    "filter",
    "fourierfilter",
    "bank",
    "morlet1d",
    "variables"
]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * ", tfile, " ...")
    include(tfile)
end
