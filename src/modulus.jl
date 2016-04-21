immutable Modulus <: AbstractPointwise
end

call(ρ::Modulus, data::AbstractArray) = abs(data)

call{N<:AbstractNode}(ρ::Modulus, pair::Pair{Path,N}) =
    (pair.first => ρ(pair.second))

call(ρ::Modulus, node::AbstractNode) = Node(ρ(node.data), node.ranges)

immutable SquaredModulus <: AbstractPointwise
end

call(ρ::SquaredModulus, data::AbstractArray) = abs2(data)
