module WaveletScattering

# using Clp
using DataStructures
using JuMP
using MathProgBase
using Mocha

include("domain.jl")
include("group.jl")
include("spec.jl")
include("meta.jl")
include("filter.jl")
include("fourier1dfilter.jl")
include("bank.jl")
include("morlet1d.jl")
include("path.jl")
include("node.jl")
include("pointwise.jl")
include("blob.jl")
include("layer.jl")
include("layerstate.jl")

end
