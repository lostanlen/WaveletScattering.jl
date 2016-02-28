# WaveletLayerState
immutable WaveletLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    blobs::Vector{B}
    blobs_diff::Vector{B}
    layer::WaveletLayer
end

function forward!(backend::Mocha.CPUBackend, state::PointwiseLayerState,
                  ρ::AbstractPointwise, inputs::Vector)
    @inbounds for idblob in eachindex(inputs)
        map!(ρ, state.blobs[idblob], inputs[idblob])
    end
end

function Mocha.setup{T,N}(
        backend::Mocha.CPUBackend,
        diffs:Vector{ScatteredBlob{T}},
        inputs::Vector{Mocha.CPUBlob{T,N}},
        layer::WaveletLayer)
    blobs = Vector{Mocha.Blob}(length(inputs))
    for idblob in eachindex(inputs)
        innodes = inputs[idblob].nodes
        outnodes = Dict{Path, AbstractFourierNode}()
        for path in keys(innodes)

        end
        blobs[idblob] = ScatteredBlob(outnodes)
    end
    blobs_diff = 0
    WaveletLayerState(blobs, blobs_diff, layer)
end
