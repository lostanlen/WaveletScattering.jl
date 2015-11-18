abstract AbstractBank{
        T<:Number,
        D<:AbstractDomain,
        G<:AbstractPointGroup,
        W<:RedundantWaveletClass}

immutable Bank1D{
        T<:Number,
        D<:LineDomains,
        G<:LineGroups,
        W<:RedundantWaveletClass,
        S<:Spec1D} <: AbstractBank{T,D,G,W}
    ϕ::AbstractFilter{T,D,G,W}
    ψs::Array{AbstractFilter{T,D,G,W},2}
    behavior::Behavior{D}
    spec::S
    function call{T,D,G,W}(
            ::Type{Bank1D}, spec::Spec1D{T,D,G,W}, behavior::Behavior{G})
        ψs = pmap(AbstractFilter, metas, fill(spec, length(metas)))
        ψs = convert(Array{AbstractFilter{T,D,G,W},2}, ψs)
        ϕ = scalingfunction(spec)
        renormalize!(ϕ, ψs, metas, spec)
        behavior = Behavior(ϕ, ψs, spec,
            is_ϕ_applied, j_range, log2_oversampling, max_log2_stride)
        new{T}(ϕ, ψs, behavior, metas, spec)
    end
end

"""A `FourierOriented1DBank` is a one-dimensional, oriented filter bank defined
in the Fourier domain. It is ""oriented"" insofar as its filters have negative
center frequencies as well as positive center frequencies, as represented by
the orientation parameter `θ`. It is advisable to use this type of filter bank
when handling complex 1d data of moderate to large length."""
immutable FourierOriented1DBank{T<:FFTW.fftwNumber} <: AbstractOrientedBank{T}
    ϕ::Symmetric1DFilter{T}
    ψs::Matrix{AbstractFourierFilter{T,1}}
    behavior::Behavior
    metas::Matrix{OrientedItem}
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
        @inbounds metas = [ OrientedItem(
            γs[γ], θs[θ], χs[γ], bws[γ], ξs[γ], js[γ], qs[γ], scs[γ])
            for γ in eachindex(γs), θ in eachindex(θs) ]
        ψs = pmap(fourierwavelet, metas[:, 1], fill(spec, length(metas)))
        ψs = convert(Array{AbstractFourierFilter{T,1}}, ψs)
        ψs = hcat(ψs, map(spin, ψs))
        ϕ = scalingfunction(spec)
        renormalize!(ϕ, ψs, metas, spec)
        behavior = Behavior(ϕ, ψs, spec,
            is_ϕ_applied, j_range, log2_oversampling, max_log2_stride)
        new{T}(ϕ, ψs, behavior, metas, spec)
    end
end
FourierOriented1DBank(spec::Abstract1DSpec ; args...) =
    FourierOriented1DBank{spec.signaltype}(spec, args...)
