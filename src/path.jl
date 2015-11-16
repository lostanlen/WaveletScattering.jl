"A `Literal` is given by a symbol and a depth > 0, defaulting to 1."
immutable Literal
    symbol::Symbol
    depth::Int
end
Literal(tup::Tuple{Symbol,Int}) = Literal(tup[1], tup[2])
Literal(sym::Symbol) = Literal(sym, 1)

"A `PathKey` is a double-ended queue (Deque) of `Literal`s"
type PathKey
    deque::DataStructures.Deque{Literal}
    PathKey() = new(DataStructures.Deque{Literal}())
    function PathKey(args...)
        deque = DataStructures.Deque{Literal}()
        idarg = 1
        nargs = length(args)
        while idarg <= nargs
            if (idarg+1)<=nargs &&
                    isa(args[idarg], Symbol) && isa(args[idarg+1], Int)
                DataStructures.push!(deque,
                    Literal(args[idarg], args[idarg+1]))
                idarg += 2
            else
                DataStructures.push!(deque, Literal(args[idarg]))
                idarg += 1
            end
        end
        new(deque)
    end
end

back(pathkey::PathKey) = DataStructures.back(pathkey.deque)

isempty(pathkey::PathKey) = DataStructures.isempty(pathkey.deque)

pop!(pathkey::PathKey) = PathKey(DataStructures.pop!(pathkey.deque))

"""A `Path` is a dictionary whose keys are `PathKey`s and whose values are
integer indices"""
immutable Path
    _dict::Dict{PathKey,Int}
    function Path(pairs...)
        _dict = Dict{PathKey,Int}()
        for pair in pairs
            push!(_dict, PathKey(pair.first...) => pair.second)
        end
        new(_dict)
    end
end

typealias PathRange Pair{PathKey,StepRange{Int,Int}}
