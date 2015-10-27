abstract AbstractScatteredLayerState <: Mocha.LayerState

# WaveletLayerState
immutable WaveletLayerState{W<:AbstractBank,B<:AbstractScatteredBlob} <:
        AbstractScatteredLayerState
    bank::W
    blobs::Vector{B}
    blobs_diff::Any
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

function setup{T<:FFTW.fftwReal,N}(
        backend::Mocha.CPUBackend,
        layer::WaveletLayer,
        bank::FourierNonOriented1DBank{T},
        inputs::Vector{Mocha.CPUBlob{T,N}},
        diffs::Vector{RealFourierBlob{T,N}} ;
        subscripts = (:time,))
    nBlobs = length(inputs)
    blobs = Array(RealFourierBlob{T,N}, nBlobs)
    subscripts = setup_subscripts(subscripts, N)
    for idblob in eachindex(inputs)
        blobs[idblob] = RealFourierBlob(inputs[idblob].data, subscripts)
    end
    blobs_diff = 0
    WaveletLayerState(bank, blobs, blobs_diff, layer)
end

function setup_subscripts{NS}(subscripts::NTuple{NS,Symbol}, ND::Int)
    suffix_subscripts = ntuple(n -> symbol(:var, n), ND-NS)
    return (subscripts..., suffix_subscripts...)
end
