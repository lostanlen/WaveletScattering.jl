module WaveletScattering

# using Clp
using DataStructures
using JuMP
using MathProgBase
using Mocha
using Wavelets

include("domain.jl")
include("group.jl")
include("waveletclass.jl")
include("path.jl")
include("meta.jl")
include("spec.jl")
include("spec1d.jl")
include("filter.jl")
include("fourier1dfilter.jl")
include("behavior.jl")
include("bank.jl")
include("morlet1d.jl")
include("node.jl")
include("blob.jl")
include("pointwise.jl")
include("layer.jl")
include("layerstate.jl")
include("fourierlayerstate.jl")
include("waveletlayerstate.jl")

end
