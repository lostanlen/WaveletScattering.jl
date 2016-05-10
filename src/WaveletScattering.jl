module WaveletScattering

# using Clp
using DataStructures
using JuMP
using MathProgBase
using Mocha

include("domain.jl")
include("group.jl")
include("waveletclass.jl")
include("path.jl")
include("meta.jl")
include("spec.jl")
include("spec1d.jl")
include("spec2d.jl")
include("filter.jl")
include("fourier1dfilter.jl")
include("weighting.jl")
include("behavior.jl")
include("morlet1d.jl")
include("bank.jl")
include("bank1d.jl")
include("bank2d.jl")
include("node.jl")
include("blob.jl")
include("pointwise.jl")
include("modulus.jl")
include("log1p.jl")
include("pointwise-call.jl")
include("symbols.jl")
include("inputlayer.jl")
include("pointwiselayer.jl")
include("fourierlayer.jl")
include("waveletlayer.jl")
include("layerstate.jl")
include("inputlayerstate.jl")
include("pointwiselayerstate.jl")
include("fourierlayerstate.jl")
include("waveletlayerstate.jl")

end
