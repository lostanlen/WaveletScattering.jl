"A `Literal` is given by a symbol and a level > 0, defaulting to 1."
immutable Literal
    symbol::Symbol
    level::Int
end
typealias SymbolInt Tuple{Symbol,Int}
Literal(tup::SymbolInt) = Literal(tup[1], tup[2])
Literal(sym::Symbol) = Literal(sym, 1)

"A `PathKey` is a linked list of `Literal`s"
type PathKey
    list::LinkedList{Literal}
    PathKey() = new(nil(Literal))
    PathKey(args...) = new(map(Literal, list(args...)))
end

"""A `Path` is a dictionary whose keys `PathKeys` and whose values are
integer indices"""
typealias Path Dict{PathKey, Int}

typealias PathRange Dict{PathKey, Range}
