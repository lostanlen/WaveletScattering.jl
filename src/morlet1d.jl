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
        new(RealT(ɛ), T, log2_length, RealT(max_qualityfactor),
            RealT(max_scale), nFilters_per_octave, nOctaves)
    end
end

function Morlet1DSpec{T<:Number}(opts::Options{CheckError})
    RealT = realtype(T)
    @defaults opts ɛ = eps(RealT)
    @defaults opts log2_length = 15
    @defaults opts max_qualityfactor = one(RealT)
    if haskey(opts.key2index, :nFilters_per_octave)
        @defaults opts max_qualityfactor = float(opts[:nFilters_per_octave])
    else
        @defaults opts max_qualityfactor = one(RealT)
    end
    @defaults opts max_scale = RealT(Inf)
    @defaults opts nFilters_per_octave = ceil(max_qualityfactor)
    @defaults opts nOctaves = 8
    @check_used opts
    checkspec(opts)
    Morlet1DSpec{precision}(opts[:ɛ], opts[:log2_length],
             opts[:max_qualityfactor], opts[:max_scale],
             opts[:nFilters_per_octave], opts[:nOctaves])
end

# Real input by default
Morlet1DSpec(opts::Options) = Morlet1DSpec(opts)
Morlet1DSpec{Number}(opts::Options) = Morlet1DSpec{Real}(opts)

# Single-precision input by default
Morlet1DSpec{Real}(opts::Options) = Morlet1DSpec{Float32}(opts)
Morlet1DSpec{Complex}(opts::Options) = Morlet1DSpec(Complex{Float32})(opts)

# Empty options by default
Morlet1DSpec() = Morlet1DSpec(@options)
Morlet1DSpec{T<:Number}() = Morlet1DSpec{T}(@options)
