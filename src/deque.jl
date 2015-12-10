function Base.(:(==))(x::Deque, y::Deque ; front::Bool = true)
    itx = DataStructures.start(x)
    ity = DataStructures.start(y)
    while ~DataStructures.done(x, itx) && ~DataStructures.done(y, ity)
        (valx, itx) = DataStructures.next(x, itx)
        (valy, ity) = DataStructures.next(y, ity)
        (valx == valy) || return false
    end
    return DataStructures.done(x, itx) && DataStructures.done(y, ity)
end
