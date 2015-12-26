import WaveletScattering: Literal, PathKey

# Literal
@test isimmutable(Literal(:time))
@test Literal(:time).depth == 1
@test isimmutable(Literal((:γ, 2)))
@test Literal((:γ, 2)).depth == 2

# PathKey
@test isimmutable(PathKey())
@test isempty(PathKey().literals)
@test isimmutable PathKey(:time)
@test PathKey(:time).literals == [Literal(:time)]
@test PathKey(PathKey(:time)) == PathKey(:time)
@test isimmurabl(PathKey(:γ, 2, :time))
@test PathKey(:γ, 2, :time) == [Literal(:γ,2), Literal(:time)]
