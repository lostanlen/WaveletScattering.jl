"""A `Meta` object contains all the meta-information to identify a wavelet
within a filter bank."""
abstract AbstractMeta

"""A `NonOrientedMeta` object contains all the meta-information to identify a
non-oriented wavelet within a filter bank. Fields:
* γ log-scale. 2^(-γ) is proportional to center frequency
* θ orientation. θ is proportional to angle (2d) or spin (1d)
* χ chroma. χ is mod(γ, nFilters_per_octave)
* j octave. j is div(γ, nFilters_per_octave)
* bandwidth ∈]0,1] is the width at -3dB, expressed in fraction of signal length
* centerfrequency ∈]0,1] is expressed in fraction of signal length
* qualityfactor ∈[1,max_qualityfactor] is equal to centerfrequency/bandwidth
* scale is the FWTM (full width at tenth maximum) in spatial domain"""
immutable NonOrientedMeta <: AbstractMeta
    γ::Int16
    χ::Int8
    bandwidth::Float64
    centerfrequency::Float64
    j::Int8
    qualityfactor::Float64
    scale::Float64
end
