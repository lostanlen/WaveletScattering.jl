immutable Morlet1DSpec{T<:Number} <: Abstract1DSpec{T}
    ɛ # ::realtype(T) (imposed by inner constructor)
    log2_length::Int
    max_qualityfactor::Float64
    max_scale::Float64
    nFilters_per_octave::Int
    nOctaves::Int
    signaltype::Type{T}
    function Morlet1DSpec{T}(signaltype::Type{T},
      ɛ=eps(realtype(T)),
      log2_length=15,
      max_qualityfactor=realtype(T)(nFilters_per_octave),
      max_scale=Inf,
      nFilters_per_octave=1,
      nOctaves=Inf)
        if isinf(nOctaves)
            log2_nFilters_per_octave = ceil(Int, log2(nFilters_per_octave))
            log2_max_qualityfactor = ceil(Int, log2(max_qualityfactor))
            if max_scale < (exp2(log2_length)+eps(RealT))
                gap = max(1+log2_nFilters_per_octave, 2)
            else
                gap = max(1+log2_nFilters_per_octave, 2+log2_max_qualityfactor)
            end
            nOctaves = log2_length - gap
        end
        checkspec(ɛ, log2_length, max_qualityfactor,
          max_scale, nFilters_per_octave, nOctaves)
        new(realtype(T)(ɛ), log2_length, max_qualityfactor,
          max_scale, nFilters_per_octave, nOctaves, signaltype)
    end
end

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
