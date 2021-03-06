using Base.Test

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
    "behavior",
    "morlet1d",
    "bank",
    "bank1d",
    "blob",
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
