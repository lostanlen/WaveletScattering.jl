using Base.Test

# path.jl
import WaveletScattering: Literal, PathKey

# Literal
@test isimmutable(Literal(:time))
@test Literal(:time).depth == 1
@test isimmutable(Literal((:γ, 2)))
@test Literal((:γ, 2)).depth == 2

# PathKey
@test isimmutable(PathKey())
@test isempty(PathKey().literals)
@test isimmutable(PathKey(:time))
@test PathKey(:time).literals == [Literal(:time)]
@test PathKey(PathKey(:time)) == PathKey(:time)
@test isimmutable(PathKey(:γ, 2, :time))
@test PathKey(:γ, 2, :time).literals == [Literal(:γ,2), Literal(:time)]
