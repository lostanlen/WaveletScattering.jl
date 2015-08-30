"A `Literal` is given by a symbol and a level > 0, defaulting to 1."
immutable Literal
    symbol::Symbol
    level::Int
end
typealias SymbolInt Tuple{Symbol,Int}
Literal(tup::SymbolInt) = Literal(tup[1], tup[2])
Literal(sym::Symbol) = Literal(sym, 1)

"A `Variable` is a linked list of `Literal`s"
typealias VariableKey LinkedList{Literal}
variablekey() = nil(Literal)
variablekey(head, tail...) = cons(Literal(head), variablekey(tail...))

"""A `VariableTree` is a hybrid `Dict`-of-`Vector`s recursive container
indexed by the scattering symbols and levels that make up a `Literal`."""
type VariableTree{V}
    value::V
    symbols::Dict{Symbol,Vector{VariableTree{V}}}
end
VariableTree{V}(value::V) =
    VariableTree(value, Dict{Symbol, Vector{VariableTree{V}}}())

"A `VariableTree` leaf can be read with a `VariableKey` accessor"
Base.getindex(t::VariableTree, key::Nil) = t.value
Base.getindex(t::VariableTree, key::Cons) =
    getindex(t.symbols[key.head.symbol][key.head.level], key.tail)

"A `VariableTree` leaf can be written with a `VariableKey` accessor"
function Base.setindex!(t::VariableTree, value, key::Nil)
    t.value = value
    return t
end
function Base.setindex!(t::VariableTree, value, key::Cons)
  setindex!(t.symbols[key.head.symbol][key.head.level], value, key.tail)
  return t
end

"Returns `true` if a `VariableTree` contains a given key, `false` otherwise."
Base.haskey(t::VariableTree, key::Nil) = true
function Base.haskey(t::VariableTree, key::Cons)
    haskey(t.symbols, key.head.symbol) || return false
    length(t.symbols[key.head.symbol])>=key.head.level || return false
    return haskey(t.symbols[key.head.symbol][key.head.level], key.tail)
end

"Traverses a `VariableTree` up to the given key, yields the subsequent branch."
subtree(t::VariableTree, key::Nil) = t
subtree(t::VariableTree, key::Cons) =
    subtree(t.symbols[key.head.symbol][key.head.level], key.tail)
