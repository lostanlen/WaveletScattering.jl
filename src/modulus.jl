immutable Modulus <: AbstractPointwise
end

call(ρ::Modulus, data::AbstractArray) = abs(data)

immutable SquaredModulus <: AbstractPointwise
end

call(ρ::SquaredModulus, data::AbstractArray) = abs2(data)
