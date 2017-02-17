immutable PointwiseLayerState{P<:AbstractPointwise} <:
        AbstractScatteredLayerState
    blobs::Vector{Mocha.Blob}
    layer::PointwiseLayer{P}
end

function PointwiseLayerState(
        backend::Mocha.CPUBackend,
        layer::PointwiseLayer,
        inputs::Vector{Mocha.Blob})
    blobs = Vector{Mocha.Blob}(length(inputs))
    for idblob in eachindex(inputs)
        blobs[idblob] = layer.ρ(inputs[idblob])
    end
    return PointwiseLayerState(blobs, layer)
end

function Mocha.setup(
        backend::Mocha.CPUBackend,
        layer::PointwiseLayer,
        inputs::Vector{Mocha.Blob},
        diffs::Vector{Mocha.Blob})
    return PointwiseLayerState(backend, layer, inputs)
end

function forward(
        backend::Mocha.CPUBackend,
        layerstate::PointwiseLayerState,
        inputs::Vector{Mocha.Blob})
    for id in eachindex(inputs)
        map!(layerstate.ρ, layerstate.blobs[id], inputs[id])
    end
end

function forward!(
        backend::Mocha.CPUBackend,
        state::PointwiseLayerState,
        ρ::AbstractPointwise,
        inputs::Vector)
    @inbounds for idblob in eachindex(inputs)
        map!(ρ, state.blobs[idblob], inputs[idblob])
    end
end
