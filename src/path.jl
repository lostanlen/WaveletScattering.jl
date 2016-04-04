"A `Literal` is given by a symbol and a depth > 0, defaulting to 1."
immutable Literal
    symbol::Symbol
    depth::Int
end
Literal(tup::Tuple{Symbol,Int}) = Literal(tup[1], tup[2])
Literal(sym::Symbol) = Literal(sym, 1)

function Base.string(literal::Literal)
    if literal.depth == 1
        return string(literal.symbol)
    else
        return string(literal.symbol) * string(literal.depth)
    end
end

Base.isless(a::Literal, b::Literal) =
    (a.symbol<b.symbol) || ((a.symbol==b.symbol) && (a.depth<b.depth))


"A `PathKey` is a `Vector` of `Literal`s."
immutable PathKey
    literals::Vector{Literal}
    PathKey() = new(Literal[])
    PathKey(pathkey::PathKey) = pathkey
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

Base.(:(==))(x::PathKey, y::PathKey) = (x.literals == y.literals)

Base.convert(::Type{PathKey}, sym::Symbol) = PathKey(sym)
Base.convert(::Type{PathKey}, tup::Tuple) = PathKey(tup...)

function Base.isless(a::PathKey, b::PathKey)
    for (x, y) in zip(a.literals, b.literals)
        if (x.symbol > y.symbol) || (x.depth > y.depth)
            return false
        end
    end
    if length(a.literals) > length(b.literals)
        return false
    end
    return a.literals != b.literals
end

function Base.string(pathkey::PathKey)
    return join([string(literal) for literal in pathkey.literals], "_")
end

Base.symbol(pathkey::PathKey) = symbol(string(pathkey))


"""A `Path` is a sorted dictionary whose keys are `PathKey`s and whose
values are integer indices."""
immutable Path
    sdict::DataStructures.SortedDict{PathKey,Int}
    function Path(pairs...)
        sdict = DataStructures.SortedDict(Pair{PathKey,Int}[])
        for pair in pairs
            push!(sdict, PathKey(pair.first) => pair.second)
        end
        new(sdict)
    end
end

Base.(:(==))(x::Path, y::Path) = (x.sdict == y.sdict)


"""A `PathRange` is a sorted dictionary whose keys are `PathKey`'s
and whose values are integer ranges."""
immutable PathRange
    sdict::DataStructures.SortedDict{PathKey,StepRange{Int,Int}}
    function PathRange(pairs...)
        sdict = DataStructures.SortedDict(Pair{PathKey,StepRange{Int,Int}}[])
        for pair in pairs
            rhs = pair.second
            if isa(rhs, Int)
                push!(sdict, PathKey(pair.first) => rhs:1:rhs)
            elseif isa(pair.second, UnitRange{Int})
                push!(sdict, PathKey(pair.first) => rhs.start:1:rhs.stop)
            elseif isa(pair.second, StepRange{Int,Int})
                push!(sdict, PathKey(pair.first) => rhs)
            end
        end
        new(sdict)
    end
end

Base.(:(==))(x::PathRange, y::PathRange) = (x.sdict == y.sdict)
