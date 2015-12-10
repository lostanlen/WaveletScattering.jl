function Base.isequal(x::Deque, y::Deque)
    DataStructures.isempty(x) && DataStructures.isempty(y) && return true
    isequal(DataStructures.pop!(x), DataStructures.pop!(y)) || return false
    isequal(x,y)
end
