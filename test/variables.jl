import WaveletScattering: Literal, variablekey, VariableTree,
                          getindex, setindex!, haskey, subtree

# Literal
time_literal = Literal(:time)
gamma2_literal = Literal((:γ, 2))
@test time_literal.level == 1
@test gamma2_literal.level == 2
@test isimmutable(time_literal)
@test isimmutable(gamma2_literal)

# variablekey
@test isa(variablekey(),Nil{Literal})
varkey = variablekey(:time, (:γ, 2))
@test varkey.head == time_literal
@test varkey.tail.head == gamma2_literal
@test varkey.tail.tail == variablekey()

# VariableTree constructor
value = 1.0
variabletree = VariableTree(value)
@test variabletree.value == value
@test isa(variabletree.symbols, Dict{Symbol,Vector{VariableTree{Float64}}})
@test isempty(variabletree.symbols)
variabletree.symbols[:time] = [VariableTree(2.0)]

# VariableTree getindex
@test variabletree[nil()] == 1.0
@test variabletree[variablekey(:time)] == 2.0

# VariableTree setindex!
variabletree[variablekey(:time)] = 3.0
@test variabletree[variablekey(:time)] == 3.0

# VariableTree haskey
@test haskey(variabletree, nil())
@test haskey(variabletree, variablekey(:time))
@test !haskey(variabletree, variablekey(:space))
@test !haskey(variabletree, variablekey((:time,2)))

# subtree
@test subtree(variabletree, nil())[nil()] == 1.0
@test subtree(variabletree, nil())[variablekey(:time)] == 3.0
@test subtree(variabletree, variablekey(:time))[nil()] == 3.0
