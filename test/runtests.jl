using WaveletScattering

tests = [
    "variables"
]

println("Running tests:")

for t in tests
    println(" * ", t, ".jl ...")
    include(t)
end
