abstract AbstractPointwise

immutable Identity <: AbstractPointwise
end

call{T,N}(ρ::Identity, data::AbstractArray{T,N}) = data

immutable Log1P{T<:AbstractFloat} <: AbstractPointwise
    threshold::T
    Log1P(threshold) = (threshold>=0 || throw(DomainError)) && new(threshold)
end

call{T,N}(ρ::Log1P, data::AbstractArray{T,N}) = log1p(ρ.threshold * data)

function Base.map(ρ::AbstractPointwise, blob_in::ScatteredBlob)
    blob_out = map(pair -> (pair.first => ρ(pair.second), blob_in))
end

function Base.map!(
        ρ::AbstractPointwise,
        blob_out::ScatteredBlob,
        blob_in::ScatteredBlob)
    @inbounds for id in eachindex(blob_in.nodes)
        map!(ρ, blob_out.nodes[id].data, blob_in.nodes[id].data)
    end
end

function Base.map!{T<:Real}(
        ρ::Log1P{T},
        data_out::Array{T},
        data_in::Array{T})
    @inbounds @fastmath for id in eachindex(data_in)
        data_out[id] = data_in[id] * ρ.threshold
        data_out[id] = log1p(data_out[id])
    end
end

for T in subtypes(AbstractPointwise)
    @eval begin
        function call{NODE<:AbstractNode}(
                ρ::T,
                pair::Pair{Path,Node})
            return pair.first => Node(ρ(pair.second.data), pair.second.ranges)
        end
    end
end
