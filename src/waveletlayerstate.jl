# WaveletLayerState
immutable WaveletLayerState{B<:ScatteredBlob} <: AbstractScatteredLayerState
    layer::WaveletLayer
    blobs::Vector{B}
    blobs_diff::Vector{B}
end

function forward!(backend::Mocha.CPUBackend, state::PointwiseLayerState,
                  ρ::AbstractPointwise, inputs::Vector)
    @inbounds for idblob in eachindex(inputs)
        map!(ρ, state.blobs[idblob], inputs[idblob])
    end
end

function Mocha.setup{
    T<:FFTW.fftwReal,
    D<:FourierDomain,
    G<:LineGroups,
    W<:RedundantWaveletClass,
    N}(
        backend::Mocha.CPUBackend,
        layer::WaveletLayer,
        bank::Bank1D{T,D,G,W},
        inputs::Vector{Mocha.CPUBlob{T,N}},
        diffs::Vector{ScatteredBlob{T}} ;
        symbols::Vector{Symbol} = [:time],
        flags = FFTW.ESTIMATE,
        timelimit = Inf)
    blobs = Array(ScatteredBlob{RealFourierNode{T,N}}, length(inputs))
    appendsymbols!(symbols, N)
    subscripts = map(PathKey, symbols)
    for idblob in eachindex(inputs)
        subscripts = inputs[idblob].subscripts
        region = find(subscripts == bank.behavior.pathkey)
        node = RealFourierNode(inputs[idblob].data, region, subscripts)
        nodes = Dict(Path() => node)
        blobs[idblob] = ScatteredBlob(nodes, symbols)
    end
    blobs_diff = 0
    WaveletLayerState(bank, blobs, blobs_diff, layer)
end
