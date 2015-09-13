"A `Literal` is given by a symbol and a level > 0, defaulting to 1."
immutable Literal
    symbol::Symbol
    level::Int
end
typealias SymbolInt Tuple{Symbol,Int}
Literal(tup::SymbolInt) = Literal(tup[1], tup[2])
Literal(sym::Symbol) = Literal(sym, 1)

"A `PathKey` is a linked list of `Literal`s"
typealias PathKey LinkedList{Literal}
PathKey() = nil(Literal)
PathKey(head, tail...) = cons(Literal(head), PathKey(tail...))
