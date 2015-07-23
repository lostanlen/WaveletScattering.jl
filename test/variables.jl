# Literal
time_literal = Literal(:time)
gamma2_literal = Literal((:γ, 2))
@test time_literal.level == 1
@test gamma2_literal.level == 2
@test isimmutable(time_literal)
@test isimmutable(gamma2_literal)

# VariableKey
@test isa(VariableKey(),Nil{Literal})
variablekey = VariableKey(:time, (:γ, 2))
@test variablekey.head == time_literal
@test variablekey.tail.head == gamma2_literal
@test variablekey.tail.tail == VariableKey()

# VariableTree constructor
value = 1.0
variabletree = VariableTree(value)
@test variabletree.value == value
@test isa(variabletree.symbols, Dict{Symbol,Vector{VariableTree{Float64}}})
@test isempty(variabletree.symbols)
variabletree.symbols[:time] = [VariableTree(2.0)]

# VariableTree getindex
@test variabletree[nil()] == 1.0
@test variabletree[VariableKey(:time)] == 2.0

# VariableTree setindex!
variabletree[VariableKey(:time)] = 3.0
@test variabletree[VariableKey(:time)] == 3.0

# VariableTree haskey
@test haskey(variabletree, nil())
@test haskey(variabletree, VariableKey(:time))
@test !haskey(variabletree, VariableKey(:space))
@test !haskey(variabletree, VariableKey((:time,2)))
