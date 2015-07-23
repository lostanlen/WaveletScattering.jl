using WaveletScattering

tests = [
    "variables.jl"
]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * ", t, ".jl ...")
    include(tfile)
end