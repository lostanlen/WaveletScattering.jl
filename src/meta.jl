abstract AbstractMeta
abstract AbstractΨMeta <: AbstractMeta

"""A `ψMeta` object contains all the meta-information to identify
an oriented wavelet within a filter bank. Fields:
* `γ` log-scale. `2^(-γ)` is proportional to center frequency
* `θ` orientation, i.e. angle (in 2d) or sign of center frequency (in 1d)
* `χ` chroma. `χ` is `mod(γ, nFilters_per_octave)`
* `j` octave. `j` is `div(γ, nFilters_per_octave)`
* `bandwidth` ∈]0,1] the width at -3dB, expressed in fraction of signal length
* `centerfrequency` ∈]0,1] is expressed in fraction of signal length
* `qualityfactor` ∈[1,max_qualityfactor] is equal to `centerfrequency/bandwidth`
* `scale` is the FWTM (full width at tenth maximum) in spatial domain"""
immutable ΨMeta1D <: AbstractΨMeta
    γ::Int16
    θ::Int8
    χ::Int8
    bandwidth::Float64
    centerfrequency::Float64
    j::Int8
    qualityfactor::Float64
    scale::Float64
end

"""A `ΦMeta` object contains all the meta-information to identify a low-pass
filter within a wavelet filter bank. Fields:
* `bandwidth` ∈]0,1] the width at -3dB, expressed in fraction of signal length
* `scale` ∈]0,1] is the FWTM (full width at tenth maximum) in spatial domain."""
immutable ΦMeta <: AbstractMeta
    bandwidth::Float64
    scale::Float64
end

get_bandwidth(meta::AbstractMeta) = meta.bandwidth
get_scale(meta::AbstractMeta) = meta.scale

get_γ(meta::AbstractΨMeta) = meta.γ
get_θ(meta::AbstractΨMeta) = meta.θ
get_χ(meta::AbstractΨMeta) = meta.χ
get_centerfrequency(meta::AbstractΨMeta) = meta.centerfrequency
get_j(meta::AbstractΨMeta) = meta.j
get_qualityfactor(meta::AbstractΨMeta) = meta.qualityfactor
