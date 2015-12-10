"A `Literal` is given by a symbol and a depth > 0, defaulting to 1."
immutable Literal
    symbol::Symbol
    depth::Int
end
Literal(tup::Tuple{Symbol,Int}) = Literal(tup[1], tup[2])
Literal(sym::Symbol) = Literal(sym, 1)

"A `PathKey` is a Vector of `Literal`s."
immutable PathKey
    literals::Vector{Literal}
    PathKey() = new(Literal[])
    PathKey(p::PathKey) = p
    function PathKey(args...)
        literals = Literal[]
        idarg = 1
        nargs = length(args)
        while idarg <= nargs
            if (idarg+1)<=nargs &&
                    isa(args[idarg], Symbol) && isa(args[idarg+1], Int)
                push!(literals, Literal(args[idarg], args[idarg+1]))
                idarg += 2
            else
                push!(literals, Literal(args[idarg]))
                idarg += 1
            end
        end
        new(literals)
    end
end

convert(PathKey, sym::Symbol) = PathKey(sym)
convert(PathKey, tup::Tuple) = PathKey(tup...)

"""A `Path` is a dictionary whose keys are `PathKey`s and whose values are
integer indices."""
immutable Path
    dict::Dict{PathKey,Int}
    function Path(pairs...)
        dict = Dict{PathKey,Int}()
        for pair in pairs
            push!(dict, convert(PathKey, pair.first) => pair.second)
        end
        new(dict)
    end
end

Base.(:(==))(x::Path, y::Path) = isequal(x, y)
Base.isequal(x::Path, y::Path) = (x.dict == y.dict)

"""A `PathRange` is a dictionary whose keys are `PathKey`s are whose values are
integer ranges."""
immutable PathRange
    _dict::Dict{PathKey,StepRange{Int,Int}}
    function PathRange(pairs...)
        _dict = Dict{PathKey,StepRange{Int,Int}}()
        for pair in pairs
            rhs = pair.second
            if isa(rhs, Int)
                push!(_dict, PathKey(pair.first...) => rhs:1:rhs)
            elseif isa(pair.second, UnitRange{Int})
                push!(_dict, PathKey(pair.first...) => rhs.start:1:rhs.stop)
            elseif isa(pair.second, StepRange{Int,Int})
                push!(_dict, PathKey(pair.first...) => rhs)
            end
        end
        new(_dict)
    end
end
