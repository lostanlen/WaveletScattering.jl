"""An `AbstractBank` is a wavelet filter bank. Filter banks are of two kinds:
* `AbstractNonOrientedBank`: no orientation variable `θ`, only a scale
variable `γ`. One-dimensional filter banks with real inputs are non-oriented.
* `AbstractOrientedBank`: wavelets have an orientation variable `θ` and a scale
variable `γ`. Two-dimensional filter banks are oriented. One-dimensional filter
banks with complex outputs are also oriented, in the sense that the sign
of the center frequency is regarded as an orientation variable with two
possible values."""

"""A `FourierNonOriented1DBank` is a one-dimensional, non-oriented filter bank
defined in the Fourier domain. It is ""non-oriented"" in the sense that only the
positive frequencies are guaranteed to be covered. Indeed, assuming that the
input signal will be real — i.e. its Fourier will be symmetric — it can be
recovered from the ""positive half"" of its Fourier spectrum. It is
advisable to use this type of filter bank when handling real 1d data of
moderate to large length."""
immutable Bank{T<:FFTW.fftwNumber,D<:AbstractDomain,G<:AbstractPointGroup}
    ϕ::AbstractFilter{T,D,G}
    ψs::Array{AbstractFilter{T,D,G},3}
    behavior::Behavior{D}
    items::Vector{AbstractItem{G}}
    spec::AbstractSpec{T,G}
    function call{T<:FFTW.fftwNumber}(::Type{FourierNonOriented1DBank{T}},
        spec::AbstractSpec{T,D}, behavior::Behavior{D,G})
        T == spec.signaltype || error("""Type parameter of
        FourierNonOriented1DBankmust must be equal to spec.signaltype""")
        γs, χs, js = gammas(spec), chromas(spec), octaves(spec)
        ξs, qs = centerfrequencies(spec), qualityfactors(spec)
        scs, bws = scales(spec), bandwidths(spec)
        @inbounds items = [ NonOrientedItem(
            γs[1+γ], χs[1+γ], bws[1+γ], ξs[1+γ], js[1+γ], qs[1+γ], scs[1+γ])
            for γ in γs ]
        ψs = pmap(fourierwavelet, items, fill(spec, length(items)))
        ψs = convert(Array{AbstractFourierFilter{T,1}}, ψs)
        ϕ = scalingfunction(spec)
        renormalize!(ϕ, ψs, items, spec)
        behavior = Behavior(ϕ, ψs, spec,
            is_ϕ_applied, j_range, log2_oversampling, max_log2_stride)
        new{T}(ϕ, ψs, behavior, items, spec)
    end
end
FourierNonOriented1DBank(spec::Abstract1DSpec ; args...) =
    FourierNonOriented1DBank{spec.signaltype}(spec ; args)

"""A `FourierOriented1DBank` is a one-dimensional, oriented filter bank defined
in the Fourier domain. It is ""oriented"" insofar as its filters have negative
center frequencies as well as positive center frequencies, as represented by
the orientation parameter `θ`. It is advisable to use this type of filter bank
when handling complex 1d data of moderate to large length."""
immutable FourierOriented1DBank{T<:FFTW.fftwNumber} <: AbstractOrientedBank{T}
    ϕ::Symmetric1DFilter{T}
    ψs::Matrix{AbstractFourierFilter{T,1}}
    behavior::Behavior
    items::Matrix{OrientedItem}
    spec::Abstract1DSpec{T}
    function call{T<:FFTW.fftwNumber}(
            ::Type{FourierOriented1DBank{T}}, spec::Abstract1DSpec ;
            is_ϕ_applied::Bool = false,
            j_range::UnitRange{Int} = 0:(spec.nOctaves-1),
            log2_oversampling::Int = 0, max_log2_stride::Int = spec.nOctaves-1)
        T == spec.signaltype || error("""Type parameter of
        FourierNonOriented1DBankmust be equal to spec.signaltype""")
        γs, χs, js = gammas(spec), chromas(spec), octaves(spec)
        ξs, qs = centerfrequencies(spec), qualityfactors(spec)
        scs, bws = scales(spec), bandwidths(spec)
        θs = 0:1
        @inbounds items = [ OrientedItem(
            γs[γ], θs[θ], χs[γ], bws[γ], ξs[γ], js[γ], qs[γ], scs[γ])
            for γ in eachindex(γs), θ in eachindex(θs) ]
        ψs = pmap(fourierwavelet, items[:, 1], fill(spec, length(items)))
        ψs = convert(Array{AbstractFourierFilter{T,1}}, ψs)
        ψs = hcat(ψs, map(spin, ψs))
        ϕ = scalingfunction(spec)
        renormalize!(ϕ, ψs, items, spec)
        behavior = Behavior(ϕ, ψs, spec,
            is_ϕ_applied, j_range, log2_oversampling, max_log2_stride)
        new{T}(ϕ, ψs, behavior, items, spec)
    end
end
FourierOriented1DBank(spec::Abstract1DSpec ; args...) =
    FourierOriented1DBank{spec.signaltype}(spec, args...)
