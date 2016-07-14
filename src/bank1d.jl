"""A `Bank1D` is a one-dimensional wavelet filter bank, parametrized by
* `T`: numeric type of input, e.g. Float32, Float64.
* `D`: transform domain. Either `FourierDomain{1}` or `SpatialDomain{1}`.
* `G`: point group. Either `TrivialGroup` or `ReflectionGroup`.
* `W`: wavelet class, e.g. `Morlet` or `Gammatone`.
Its fields are
* `ϕ`: low-pass filter, also called scaling function.
* `ψs`: 3d array of wavelets, indexed by spin, chroma, and octave.
* `behavior`: mutable behavior, e.g. oversampling.
* `spec`: immutable specifications, e.g. number of filters per octave.
To create a `Bank1D`
1. define a `Spec1D`,
2. define a `PathKey` on which to apply the wavelet filterbank,
3. if needed, provide behavior characteristics as keyword arguments.
Example:
spec = Spec1D(nFilters_per_octave = 12, nOctaves = 10)
pathkey = PathKey(:time)
bank = Bank1D(spec, pathkey, j_range = 2:9)"""
immutable Bank1D{
        T<:Number,
        D<:LineDomains,
        G<:LineGroups,
        W<:RedundantWaveletClass} <: AbstractBank{T,D,G,W}
    ϕ::AbstractFilter{T,D}
    ψs::Array{AbstractFilter{T,D},3}
    behavior::Behavior{T}
    spec::Spec1D{T,D,G,W}
    function (::Type{Bank1D}){T,D,G,W}(
            spec::Spec1D{T,D,G,W},
            pathkey::PathKey = PathKey(:time) ;
            is_ϕ_applied::Bool = false,
            j_range::UnitRange{Int} = 0:(spec.nOctaves-1),
            log2_oversampling::Int = 0,
            max_log2_stride::Int = spec.nOctaves-1,
            weighting::AbstractWeighting = EqualWeighting())
        (nΘs, nΧs, nJs) = size(spec.ψmetas)
        ψs = Array(AbstractFilter{T,D}, (nΘs, nΧs, nJs))
        ψs[1, :, :] =
            pmap(AbstractFilter, spec.ψmetas[1, :, :], fill(spec, nΧs * nJs))
        (nΘs > 1) && spin!(ψs)
        ϕ = AbstractFilter(spec.ϕmeta, spec)
        renormalize!(ϕ, ψs, spec)
        behavior = Behavior(ϕ, ψs, spec, is_ϕ_applied, j_range,
            log2_oversampling, max_log2_stride, pathkey, weighting)
        new{T,D,G,W}(ϕ, ψs, behavior, spec)
    end
end

function call{T<:Real,DIM}(
        bank::Bank1D{T,FourierDomain{1},TrivialGroup},
        x::AbstractArray{T,DIM};
        flags = FFTW.ESTIMATE,
        timelimit = Inf,
        verbose = 0)
    syms = appendsymbols(fill(symbol(bank.behavior.pathkey), 1), DIM)
    inputnode = Node(x, ntuple(k -> kthrange(syms, x, k), DIM))
    fouriernode = AbstractFourierNode(inputnode, [1], flags, timelimit)
    chromakey = prepend(:χ, bank.behavior.pathkey)
    chromarange = chromakey => 0:1:(bank.spec.nFilters_per_octave-1)
    waveletranges = (collect(inputnode.ranges)..., chromarange)
    waveletnodes =
        DataStructures.SortedDict(Pair{Path,Node{Complex{T},DIM+1}}[])
    octavekey = prepend(:j, bank.behavior.pathkey)
    for j in bank.behavior.j_range
        ψ_log2_sampling = bank.behavior.ψ_log2_samplings[1+j]
        downsampled_length = size(x, 1) >> (-ψ_log2_sampling)
        octave_size = (downsampled_length,
            size(x)[2:end]..., bank.spec.nFilters_per_octave)
        octave_ft = zeros(Complex{T}, octave_size)
        inds = [fill(Colon(), DIM) ; 0]
        for χ in 0:(bank.spec.nFilters_per_octave-1)
            ψ = bank.ψs[1, 1+χ, 1+j]
            inds[end] = 1 + χ
            transform!(sub(octave_ft, inds...), ψ, fouriernode, 1)
        end
        ifft_plan = plan_ifft(octave_ft, 1 ; flags=flags, timelimit=timelimit)
        waveletnode = Node(ifft_plan * octave_ft, waveletranges)
        waveletpath = Path(octavekey => j)
        waveletnodes[waveletpath] = waveletnode
    end
    waveletblob = ScatteredBlob(waveletnodes)
    return waveletblob
end

"""Collects band-pass filters `ψ`s and low-pass filter `ϕ` in the Fourier domain. The band-pass filters are organized as a 4-D tensor, where the
dimensions correspond to frequency `ω`, orientation `θ`, chroma `χ`,
and octave index `j`."""
function Base.collect{T}(bank::Bank1D{T,FourierDomain{1}})
    N = 1 << bank.spec.log2_size
    halfN = N >> 1
    (nΘs, nΧs, nJs) = size(bank.spec.ψmetas)
    ψs = zeros(T, N, nΘs, nΧs, nJs)
    ωs = (-halfN):(halfN-1)
    for j in 0:(nJs-1), χ in 0:(nΧs-1), θ in 0:(nΘs-1)
        ψs[:, 1+θ, 1+χ, 1+j] =
            fftshift(bank.ψs[1+θ, 1+χ, 1+j][ωs])
    end
    ϕ = fftshift(bank.ϕ[ωs])
    return (ψs, ϕ)
end

function littlewoodpaleyplot{T}(bank::Bank1D{T,FourierDomain{1}})
    (ψs, ϕ) = collect(bank)
    ψs = ψs[:, :]
    ωs = linspace(0, 2pi, length(ϕ)+1)[1:(end-1)]
    lp2 = squeeze(sum(abs2(ψs), 2), 2)
    symmetrize!(lp2)
    lp2 += abs2(ϕ)
    lp = sqrt(lp2)
    Winston.plot(collect(ωs), ψs,
                 collect(ωs), ϕ,
                 collect(ωs), lp)
    Winston.xlabel("ω (rad / s)")
end

Base.ndims(::Bank1D) = 1
