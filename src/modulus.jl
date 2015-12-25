immutable Modulus <: AbstractPointwise
end

call{T,N}(ρ::Modulus, data::AbstractArray{T,N}) = abs(data)

immutable SquaredModulus <: AbstractPointwise
end

call{T,N}(ρ::SquaredModulus, data::AbstractArray{T,N}) = abs2(data)


function Base.map!{T<:Real}(
        ρ::Modulus,
        data_out::Array{T},
        data_in::Array{Complex{T}})
    map!(abs, data_out, data_in)
end
