import WaveletScattering: Literal, PathKey

# Literal
time_literal = Literal(:time)
gamma2_literal = Literal((:γ, 2))
@test time_literal.level == 1
@test gamma2_literal.level == 2
@test isimmutable(time_literal)
@test isimmutable(gamma2_literal)

# PathKey
@test isa(PathKey(),Nil{Literal})
pathkey = PathKey(:time, (:γ, 2))
@test pathkey.head == time_literal
@test pathkey.tail.head == gamma2_literal
@test pathkey.tail.tail == PathKey()
