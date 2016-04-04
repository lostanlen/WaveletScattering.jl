using Base.Test

# path.jl
import WaveletScattering: Literal, PathKey, Path

# Literal
@test isimmutable(Literal(:time))
@test Literal(:time).depth == 1
@test isimmutable(Literal((:γ, 2)))
@test Literal((:γ, 2)).depth == 2

# Conversion to String
@test string(Literal(:time)) == "time"
@test string(Literal(:γ, 2)) == "γ2"

# PathKey
@test isimmutable(PathKey())
@test isempty(PathKey().literals)
@test isimmutable(PathKey(:time))
@test PathKey(:time).literals == [Literal(:time)]
@test PathKey(PathKey(:time)) == PathKey(:time)
@test isimmutable(PathKey(:γ, 2, :time))
@test PathKey(:γ, 2, :time).literals == [Literal(:γ,2), Literal(:time)]

# Conversion from Symbol
@test convert(PathKey, :time)== PathKey(:time)

# Conversion from Tuple
@test convert(PathKey, (:γ, 2, :time)) == PathKey(:γ, 2, :time)

# Conversion to string
@test string(PathKey(:time)) == "time"
@test string(PathKey(:γ, 1, :time)) == "γ_time"
@test string(PathKey(:γ, 2, :time)) == "γ2_time"

# Conversion to symbol
@test symbol(PathKey(:time)) == :time
@test symbol(PathKey(:γ, 1, :time)) == :γ_time
@test symbol(PathKey(:γ, 2, :time)) == :γ2_time
