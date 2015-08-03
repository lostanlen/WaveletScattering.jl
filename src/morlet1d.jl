immutable Morlet1DSpec{T<:Number} <: Abstract1DSpec{T}
    ɛ # ::realtype(T) (imposed by inner constructor)
    log2_length::Int
    max_qualityfactor::Float64
    max_scale::Float64
    nFilters_per_octave::Int
    nOctaves::Int
    signaltype::Type{T}
    function Morlet1DSpec(ɛ, log2_length::Int, max_qualityfactor, max_scale,
                          nFilters_per_octave::Int, nOctaves::Int, signaltype)
        RealT = realtype(T)
        checkspec(ɛ, log2_length, max_qualityfactor, max_scale,
            nFilters_per_octave, nOctaves)
        new(RealT(ɛ), log2_length, max_qualityfactor,
            max_scale, nFilters_per_octave, nOctaves, signaltype)
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
    if opts[:nFilters_per_octave] == nothing
        @defaults opts max_qualityfactor = one(RealT)
    else
        @defaults opts max_qualityfactor = float(opts[:nFilters_per_octave])
    end
    @defaults opts nFilters_per_octave = ceil(Int, max_qualityfactor)
    # By default, no upper limit on the scale of the lowest-frequency wavelets.
    @defaults opts max_scale = RealT(Inf)
    # By default, the filter bank covers the whole frequency range, while
    # ensuring that wavelet scales remain below 2^log2_length
    log2_nFilters_per_octave = floor(Int, log2(nFilters_per_octave))
    log2_max_qualityfactor = floor(Int, log2(max_qualityfactor))
    if max_scale < (exp2(log2_length)+eps(RealT))
        gap = 1 + log2_nFilters_per_octave
    else
        gap = max(1 + log2_nFilters_per_octave, 2 + log2_max_qualityfactor)
    end
    @defaults opts nOctaves = log2_length - gap
    @check_used opts
    Morlet1DSpec{T}(ɛ, log2_length, max_qualityfactor, max_scale,
                 nFilters_per_octave, nOctaves, signaltype)
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
    return (bandwidths, centerfrequencies, qualityfactors, scales)
end
