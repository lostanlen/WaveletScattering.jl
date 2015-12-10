Base.isequal(x::Deque, y::Deque) = isequal(pop!(x), pop!(y)) && isequal(x,y)
