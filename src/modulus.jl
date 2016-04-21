immutable Modulus <: AbstractPointwise
end

call(ρ::Modulus, data::AbstractArray) = abs(data)

call{N<:AbstractNode}(ρ::Modulus, pair::Pair{Path,N}) = ρ(pair.second)

call(ρ::Modulus, node::AbstractNode) = ρ(node.data)

immutable SquaredModulus <: AbstractPointwise
end

call(ρ::SquaredModulus, data::AbstractArray) = abs2(data)
