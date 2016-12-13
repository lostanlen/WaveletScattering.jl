immutable Modulus <: AbstractPointwise
end

(ρ::Modulus)(data::AbstractArray) = abs.(data)

(ρ::Modulus)(blob_in::ScatteredBlob) = ScatteredBlob(map(ρ, blob_in.nodes))

immutable SquaredModulus <: AbstractPointwise
end

(ρ::SquaredModulus)(data::AbstractArray) = abs2.(data)

(ρ::SquaredModulus)(blob_in::ScatteredBlob) =
    ScatteredBlob(map(ρ, blob_in.nodes))
