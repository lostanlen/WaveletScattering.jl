"A `Literal` is given by a symbol and a depth > 0, defaulting to 1."
immutable Literal
    symbol::Symbol
    depth::Int
end
typealias SymbolInt Tuple{Symbol,Int}
Literal(tup::SymbolInt) = Literal(tup[1], tup[2])
Literal(sym::Symbol) = Literal(sym, 1)

"A `PathKey` is a double-ended queue (Deque) of `Literal`s"
type PathKey
    deque::DataStructures.Deque{Literal}
    PathKey() = new(DataStructures.Deque{Literal}())
    function PathKey(args...)
        deque = DataStructures.Deque{Literal}()
        for arg in args
            DataStructures.unshift!(deque, Literal(arg))
        end
        new(deque)
    end
end

"""A `Path` is a dictionary whose keys are `PathKey`s and whose values are
integer indices"""
typealias Path Dict{PathKey, Int}

typealias PathRange Dict{PathKey, StepRange{Int, Int}}
