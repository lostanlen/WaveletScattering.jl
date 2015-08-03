immutable Morlet1DSpec{T<:Number} <: Abstract1DSpec{T}
    ɛ # ::realtype(T) (imposed by inner constructor)
    signaltype::Type{T}
    log2_length::Int
    max_qualityfactor # ::realtype(T) (imposed by inner constructor)
    max_scale # ::realtype(T) (imposed by inner constructor)
    nFilters_per_octave::Int
    nOctaves::Int
    function Morlet1DSpec(ɛ, signaltype::Type{T}, log2_length::Int,
                                     max_qualityfactor, max_scale,
                                     nFilters_per_octave::Int, nOctaves::Int)
        RealT = realtype(T)
        checkspec(max_qualityfactor, nFilters_per_octave)
        new(RealT(ɛ), T, log2_length, RealT(max_qualityfactor),
            RealT(max_scale), nFilters_per_octave, nOctaves)
    end
end

function Morlet1DSpec(opts::Options{CheckError})
    # Single-precision real input by default
    @defaults opts signaltype=Float32
    T = signaltype
    RealT = realtype(T)
    # Default threshold is equal to floating-point machine precision, i.e. the
    # distance between 1.0 and the next numbe of type RealT.
    @defaults opts ɛ = eps(RealT)
    # Default window length is 32768 samples, i.e. about 740 ms at 44.1kHz
    @defaults opts log2_length = 15
    # max_qualityfactor and nFilters_per_octave are mutual defaults.
    # If none is present, both are set to one.
    @defaults opts max_qualityfactor = one(RealT)
    if opts[:nFilters_per_octave] == nothing
        @defaults opts max_qualityfactor = one(RealT)
    else
        @defaults opts max_qualityfactor = float(opts[:nFilters_per_octave])
    end
    @defaults opts nFilters_per_octave = ceil(Int, max_qualityfactor)
    # By default, no upper limit on the scale of the lowest-frequency wavelets.
    @defaults opts max_scale = RealT(Inf)
    # By default, the filter bank covers the whole frequency range, while
    # ensuring that wavelet scales remain below 2^(log2_length-1)
    gap_a = ceil(Int, log2(nFilters_per_octave)) + 1
    scale_multiplier = sqrt(log(RealT(10.0))/log(RealT(2.0)))
    mother_centerfrequency = nFilters_per_octave==1 ? RealT(0.39) :
        RealT(inv(3.0 - exp2(-1.0/nFilters_per_octave)))
    mother_scale = scale_multiplier / mother_centerfrequency
    gap_b = 1 + Int(div(log2(mother_scale) - 1.0/nFilters_per_octave, 1))
    minimum_gap = max(gap_a, gap_b)
    @defaults opts nOctaves = log2_length - minimum_gap
    if (log2_length-nOctaves)<minimum_gap
        error("Wavelet support is too large in low frequencies.\n",
        "Either increase log2_length or decrease nOctaves.\n",
        "log2_length = ", log2_length, "\n",
        "nOctaves = ", nOctaves, "\n",
        "The gap should be at least ", minimum_gap, ".")
    end
    @check_used opts
    Morlet1DSpec{T}(ɛ, signaltype, log2_length, max_qualityfactor, max_scale,
                 nFilters_per_octave, nOctaves)
end

# A zero-argument constructor falls back to default options, see above
Morlet1DSpec() = Morlet1DSpec(@options)

function localize{T<:Number}(spec::Morlet1DSpec{T})
    RealT = realtype(T)
    mother_centerfrequency = spec.nFilters_per_octave==1 ? RealT(0.39) :
        RealT(inv(3.0 - exp2(-1.0/spec.nFilters_per_octave)))
    nΓs = spec.nFilters_per_octave * spec.nOctaves
    resolutions =
        exp2(range(zero(RealT), -one(T)/spec.nFilters_per_octave, nΓs))
    centerfrequencies = mother_centerfrequency * resolutions
    scale_multiplier = sqrt(log(RealT(10.0))/log(RealT(2.0)))
    unbounded_scales =
        scale_multiplier * (spec.max_qualityfactor./centerfrequencies)
    scales = min(unbounded_scales, spec.max_scale)
    unbounded_qualityfactors = scales .* centerfrequencies / scale_multiplier
    qualityfactors = clamp(unbounded_qualityfactors, 1.0, spec.max_qualityfactor)
    bandwidths = resolutions ./ qualityfactors
    scales = scale_multiplier * qualityfactors./centerfrequencies
    if scales[end]>exp2(spec.log2_length)
        error("""Wavelet support is too large in low frequencies.
        Either increase log2_length or decrease nOctaves""")
    end
    return (bandwidths, centerfrequencies, qualityfactors, scales)
end
