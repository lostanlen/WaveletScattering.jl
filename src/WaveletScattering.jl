module WaveletScattering

# using Clp
using DataStructures
using JuMP
using MathProgBase
using Mocha

include("spec.jl")
include("meta.jl")
include("filter.jl")
include("fourierfilter.jl")
include("bank.jl")
include("morlet1d.jl")
include("path.jl")
include("node.jl")
include("blob.jl")
include("layer.jl")

end
