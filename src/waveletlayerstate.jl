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

"""diffs and inputs must contain `ScatteredBlob`'s."""
function Mocha.setup{T,N}(
        backend::Mocha.Backend,
        layer::WaveletLayer,
        inputs::Vector{Mocha.Blob{T}},
        diffs::Vector{Mocha.Blob{T}})
    blobs = Vector{Mocha.Blob}(length(inputs))
    pathkey = layer.bank.behavior.pathkey
    for idblob in eachindex(inputs)
        innodes = inputs[idblob].nodes
        # TODO: tighten typing
        outnodes =
            SortedDict{Path,AbstractFourierNode,Base.Order.ForwardOrdering}()
        for inpath in keys(innodes)
            innode = innodes[inpath]

        end
        blobs[idblob] = ScatteredBlob(outnodes)
    end
    blobs_diff = Vector{Mocha.Blob}(length(inputs))
    WaveletLayerState(blobs, blobs_diff, layer)
end
