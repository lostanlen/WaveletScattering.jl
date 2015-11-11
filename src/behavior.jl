"""A `Behavior` object contains the mutable information in a filter bank. The
values of these fields may be changed *after* the construction of the filter
bank, without having to recompute the underlying architecture. Fields:

* `is_ϕ_applied`: true if and only if the lowpass filter `ϕ` (also known as
scaling function) in addition to the bandpass filters `ψ`'s.

* `j_range`: range of wavelet octaves that are actually used in the
convolutions. Default is all of them.

* `log2_oversampling`: base-2 logarithm of the oversampling factor with respect
to the critical sampling rate. Must be positive. Default is 0, i.e. no
oversampling.

* `max_log2_stride`: base-2 logarithm of the maximum distance between
neighboring coefficients after subsampling (also known as *hop size* or
*stride*). Must be positive. Default is the number of octaves minus one, which
imposes no oversampling per se."""
type Behavior
    ϕ_log2_sampling::Int
    ψ_log2_samplings::Vector{Int}
    is_ϕ_applied::Bool
    j_range::UnitRange
    log2_oversampling::Int
    max_log2_stride::Int
end

is_ϕ_applied::Bool = false,
j_range::UnitRange{Int} = 0:(spec.nOctaves-1),
log2_oversampling::Int = 0, max_log2_stride::Int = spec.nOctaves-1

function Behavior{D}(ϕ::Symmetric1DFilter,
        ψs::AbstractArray{AbstractFourierFilter, spec::AbstractSpec,
        domain::D;
        is_ϕ_applied::Bool = false,
        j_range::UnitRange{Int} = 0:(spec.nOctaves-1),
        log2_oversampling::Int = 0,
        max_log2_stride::Int = spec.nOctaves-1)
    ϕ_critical_log2_sampling = critical_log2_sampling(ϕ, spec)
    ϕ_log2_sampling =
        clamp(ϕ_critical_log2_sampling + log2_oversampling, -max_log2_stride, 0)
    edge_ids = (j_range+1) * spec.nFilters_per_octave - 1
    ψ_critical_log2_samplings =
        Int[ critical_log2_sampling(ψ, spec) for ψ in ψs[edge_ids, 1] ]
    ψ_log2_samplings = clamp(ψ_critical_log2_samplings + log2_oversampling,
        -max_log2_stride, 0)
    max_log2_stride = - min(ϕ_log2_sampling, minimum(ψ_log2_samplings))
    Behavior(ϕ_log2_sampling, ψ_log2_samplings, is_ϕ_applied, j_range,
        log2_oversampling, max_log2_stride)
end