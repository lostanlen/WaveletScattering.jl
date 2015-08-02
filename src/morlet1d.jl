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
    # By default, the filter bank covers the whole frequency range
    @defaults opts nOctaves = log2_length - 1
    @check_used opts
    Morlet1DSpec{T}(ɛ, signaltype, log2_length, max_qualityfactor, max_scale,
                 nFilters_per_octave, nOctaves)
end

# A zero-argument constructor falls back to default options, see above
Morlet1DSpec() = Morlet1DSpec(@options)
