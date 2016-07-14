immutable Modulus <: AbstractPointwise
end

Base.call(ρ::Modulus, data::AbstractArray) = abs(data)

(ρ::Modulus)(blob_in::ScatteredBlob) = map(ρ, ScatteredBlob(blob_in.nodes))

Base.abs(blob::ScatteredBlob) = Modulus()(blob)

immutable SquaredModulus <: AbstractPointwise
end

(ρ::SquaredModulus)(data::AbstractArray) = abs2(data)

(ρ::SquaredModulus)(blob_in::ScatteredBlob) =
    map(ρ, ScatteredBlob(blob_in.nodes))
