module WaveletScattering

using Compat
using DataStructures
using JuMP
using Regression

include("spec.jl")
include("meta.jl")
include("filter.jl")
include("fourierfilter.jl")
include("bank.jl")
include("morlet1d.jl")
include("path.jl")

end
