# A scattering Literal is given by a symbol and a level > 0.
# The default level is 1.
immutable Literal
    symbol::Symbol
    level::Int
end
typealias SymbolInt @compat(Tuple{Symbol,Int})
Literal(tup::SymbolInt) = Literal(tup[1], tup[2])
Literal(sym::Symbol) = Literal(sym, 1)

# A scattering Variable is a linked list of Literals.
typealias VariableKey LinkedList{Literal}
variablekey() = nil(Literal)
variablekey(head, tail...) = cons(Literal(head), variablekey(tail...))

# A scattering VariableTree is a hybrid Dict-of-Vectors recursive container
# indexed by the scattering symbols and levels that make up a Literal.
type VariableTree{V}
    value::V
    symbols::Dict{Symbol,Vector{VariableTree{V}}}
end
VariableTree{V}(value::V) =
    VariableTree(value, Dict{Symbol, Vector{VariableTree{V}}}())

# A VariableTree leaf can be read/written with a VariableKey accessor
getindex(t::VariableTree, key::Nil) = t.value
getindex(t::VariableTree, key::Cons) =
    getindex(t.symbols[key.head.symbol][key.head.level], key.tail)
function setindex!(t::VariableTree, value, key::Nil)
    t.value = value
    return t
end
function setindex!(t::VariableTree, value, key::Cons)
  setindex!(t.symbols[key.head.symbol][key.head.level], value, key.tail)
  return t
end

# The haskey function is extended to VariableTree with VariableKey indexing
import Base.haskey
haskey(t::VariableTree, key::Nil) = true
function haskey(t::VariableTree, key::Cons)
    haskey(t.symbols, key.head.symbol) || return false
    length(t.symbols[key.head.symbol])>=key.head.level || return false
    return haskey(t.symbols[key.head.symbol][key.head.level], key.tail)
end

# Unlike getindex, which yields a value, the function subtree yields a branch
subtree(t::VariableTree, key::Nil) = t
subtree(t::VariableTree, key::Cons) =
    subtree(t.symbols[key.head.symbol][key.head.level], key.tail)
