"""A `Behavior` object contains the mutable information in a filter bank. The
values of these fields may be changed *after* the construction of the filter
bank, without having to recompute the underlying architecture. Fields:

* `is_ϕ_applied`: true if and only if the lowpass filter `ϕ` (also known as scaling function) in addition to the bandpass filters `ψ`'s.

* `j_range`: range of wavelet octaves that are actually used in the convolutions. Default is all of them.

* `log2_oversampling`: base-2 logarithm of the oversampling factor with respect to the critical sampling rate. Must be positive. Default is 0, i.e. no oversampling.

* `max_log2_stride`: base-2 logarithm of the maximum distance between neighboring coefficients after subsampling (also known as *hop size* or *stride*). Must be positive. Default is the number of octaves minus one, which imposes no oversampling per se.

* `pathkey`: key of the variable over which the filter bank is applied.

* `weighting`: scalar weighting of each wavelet frequency band. Default is `EqualWeighting()`, i.e. all weights are one. Can alternatively be set to `LoudnessWeighting(samplerate)` to model the relative loudness perceived by the human ear, as defined by the international standard 61672:2003.

* `weights`: vector of floating-point weights corresponding to the frequencies of `ψ`'s. The weight assigned to `ϕ` is the same as the weight of the `ψ` with the lowest center frequency."""
type Behavior{T<:Real}
    ϕ_log2_sampling::Int
    ψ_log2_samplings::Vector{Int}
    is_ϕ_applied::Bool
    j_range::UnitRange
    log2_oversampling::Int
    max_log2_stride::Int
    pathkey::PathKey
    weighting::AbstractWeighting
    weights::Array{T, 3}
end

"""Given a lowpass filter `ϕ`, an array of wavelets `ψs`, and the corresponding
filter bank specification `spec`, the `Behavior` outer constructor computes
the critical sampling rates of all `ψs` and `ϕ`. Then, with the inputs
`log2_oversampling` and `max_log2_stride`, these sampling rates are bounded:
* from below, by the inverse of the maximum stride.
* from above, by 1, to avoid unnecessary upsampling with respect to the original
  sample rate.
NB: if `log2_oversampling` and `max_log2_stride` are left as defaults, all
sample rates are critical."""
function Behavior{
    T<:Number,
    D<:AbstractDomain,
    G<:AbstractPointGroup,
    W<:RedundantWaveletClass}(
        ϕ::AbstractFilter{T,D},
        ψs::AbstractArray{AbstractFilter{T,D},3},
        spec::AbstractSpec{T,D,G,W},
        is_ϕ_applied::Bool,
        j_range::UnitRange{Int},
        log2_oversampling::Int,
        max_log2_stride::Int,
        pathkey::PathKey,
        weighting::AbstractWeighting)
    ϕ_critical_log2_sampling = critical_log2_sampling(ϕ, spec)
    ϕ_log2_sampling =
        clamp(ϕ_critical_log2_sampling + log2_oversampling, -max_log2_stride, 0)
    ψ_critical_log2_samplings =
        Int[ critical_log2_sampling(ψ, spec) for ψ in ψs[1, end, :] ]
    ψ_log2_samplings = clamp(ψ_critical_log2_samplings + log2_oversampling,
        -max_log2_stride, 0)
    max_log2_stride = - min(ϕ_log2_sampling, minimum(ψ_log2_samplings))
    ξs = map(get_centerfrequency, spec.ψmetas)
    weights = map(T, weight_frequencies(weighting, ξs))
    Behavior(ϕ_log2_sampling, ψ_log2_samplings, is_ϕ_applied, j_range,
        log2_oversampling, max_log2_stride, pathkey, weighting, weights)
end
