call{T,N}(ρ::Log1P{T}, data::AbstractArray{T,N}) = log1p(ρ.threshold * data)

function Base.map!{T<:Real}(
        ρ::Log1P{T},
        data_out::Array{T},
        data_in::Array{T})
    @inbounds @fastmath for id in eachindex(data_in)
        data_out[id] = data_in[id] * ρ.threshold
        data_out[id] = log1p(data_out[id])
    end
end
