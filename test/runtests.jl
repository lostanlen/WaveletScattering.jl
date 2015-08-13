using Base.Test

using DataStructures

tests = [
    "spec",
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
