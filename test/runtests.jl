using Base.Test

using DataStructures

tests = [
    "spec",
    "variables"
]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * ", tfile, " ...")
    include(tfile)
end
