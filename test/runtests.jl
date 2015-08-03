using Base.Test

using DataStructures
using OptionsMod

tests = [
    "spec",
    "morlet1d",
    "variables"
]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * ", tfile, " ...")
    include(tfile)
end
