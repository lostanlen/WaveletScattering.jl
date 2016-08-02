using Base.Test

# path.jl
import WaveletScattering: Literal, PathKey, Path, PathRange

# Literal
@test isimmutable(Literal(:time))
@test Literal(:time).depth == 1
@test isimmutable(Literal((:γ, 2)))
@test Literal((:γ, 2)).depth == 2

# Conversion to String
@test string(Literal(:time)) == "time"
@test string(Literal(:γ, 2)) == "γ2"

# isless for Literal's
@test Literal(:a) < Literal(:b)
@test Literal(:a, 1) < Literal(:b)
@test Literal(:a, 1) < Literal(:a, 2)
@test !(Literal(:a) < Literal(:a))
@test !(Literal(:a, 2) < Literal(:a, 1))
@test !(Literal(:b) < Literal(:a))

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

# isless for PathKey's
@test PathKey(:a) < PathKey(:b)
@test PathKey(:a) < PathKey(:a, :b)
@test PathKey(:a) < PathKey((:a, 2))
@test !(PathKey(:a, :b) < PathKey(:a))

# Conversion to string
@test string(PathKey(:time)) == "time"
@test string(PathKey(:γ, 1, :time)) == "γ_time"
@test string(PathKey(:γ, 2, :time)) == "γ2_time"

# Conversion to symbol
@test Symbol(PathKey(:time)) == :time
@test Symbol(PathKey(:γ, 1, :time)) == :γ_time
@test Symbol(PathKey(:γ, 2, :time)) == :γ2_time

# Path
@test Path((:a, 2) => 2, :b => 3, :a => 1) ==
    Path(:a => 1, (:a, 2) => 2, :b => 3)

# isless for Path
@test Path(:a => 1) < Path((:a, 2) => 2)
@test !(Path((:a, 2) => 1) < Path(:a => 2))
@test Path(:a => 1) < Path(:a => 2)

# PathRange
@test PathRange((:a, 2) => 2:4, :b => 3, :a => 1:1:3) ==
    PathRange(:a => 1:3, (:a, 2) => 2:1:4, :b => 3:3)
