abstract AbstractScatteredLayerState <: Mocha.LayerState

# WaveletLayerState
immutable WaveletLayerState{W<:AbstractBank,B<:AbstractScatteredBlob} <:
        AbstractScatteredLayerState
    bank::W
    blobs::Vector{B}
    blobs_diff::Vector{B}
    layer::WaveletLayer
end

function forward!{IN<:AbstractFourierBlob,OUT<:AbstractFourierBlob}(
        backend::Mocha.CPUBackend, state::WaveletLayerState{OUT},
        inputs::Vector{IN})
    for idblob in eachindex(inputs)
        fft!(input[idblob])
        forward!(backend, state.blobs[idblob], state.bank, input[idblob])
        ifft!(state.blobs[idblob])
    end
end

function forward!(backend::Mocha.CPUBackend, state::WaveletLayerState,
                  ρ::AbstractPointwise, inputs::Vector)
    @inbounds for idblob in eachindex(inputs)
        map!(ρ, state.blobs[idblob], inputs[idblob])
    end
end
