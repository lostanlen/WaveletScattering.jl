abstract AbstractScatteredLayerState <: Mocha.LayerState

# WaveletLayerState
immutable WaveletLayerState{BANK<:AbstractBank,BLOB<:ScatteredBlob} <:
        AbstractScatteredLayerState
    bank::BANK
    blobs::Vector{BLOBS}
    blobs_diff::Any
    layer::WaveletLayer
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
        diffs::Vector{ScatteredBlob{T,N}} ;
        symbols::Vector{Symbol} = [:time], fourierdims::Tuple{Int} = (1,))
    blobs = Array(ScatteredBlob{RealFourierNode{T,N}}, length(inputs))
    appendsymbols!(symbols, N)
    subscripts = map(PathKey, symbols)
    for idblob in eachindex(inputs)
        node = AbstractFourierNode(inputs[idblob].data, fourierdims, subscripts)
        blobs[idblob] = RealFourierBlob(inputs[idblob].data, symbols)
    end
    blobs_diff = 0
    WaveletLayerState(bank, blobs, blobs_diff, layer)
end

function appendsymbols!(symbols::Vector{Symbol}, N::Int)
    symbols = vcat(symbols, [ symbol(:var, n) for n in 1:(length(symbols)-N) ])
end
