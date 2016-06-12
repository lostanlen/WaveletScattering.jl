function AbstractFilter{T<:FFTW.fftwReal,G<:LineGroups}(ψmeta::ΨMeta,
        spec::Spec1D{T,FourierDomain{1},G,Morlet})
    """1. **Gaussian denominator `den = 2σ²`**
    The FWHM (full width at half maximum) bw of a Gaussian bell curve of
    variance `σ` is defined by the equation
        `g(±bw/2) = exp(- (bw/2)²/(2σ²)) = 1/sqrt(2)`
    which leads to
        `bw² = 2 log(2) * 2σ²`.
    The denominator `den = 2σ²` of the Gaussian is thus equal to
        `den = 2σ² = bw² / (2 log(2))`."""
    @inbounds log2_length = spec.log2_size
    N = 1 << log2_length
    halfN = N >> 1
    center = N * T(ψmeta.centerfrequency)
    bw = N * T(ψmeta.bandwidth)
    den = @fastmath bw * bw / T(2.0 * log(2.0))
    """2. **Angle**"""
    nOrientations = get_nOrientations(spec.pointgroup)
    θ = T(2π * ψmeta.θ / nOrientations)
    """3. **Number of periods**"""
    halfsupport = sqrt(den * log(inv(spec.ɛ)))
    firstω = max(center - halfsupport, -5N/2)
    lastω = min(center + halfsupport, 5N/2 - 1)
    nPeriods = 1 + ceil(Int, (lastω-halfN) / N)
    #"""3. **Call to morlet**"""
    aspectratio = T(ψmeta.aspectratio)
    y = morlet(FourierDomain(1), θ, aspectratio, center, den, N, nPeriods)
    #"""4. **Trimming to true support boundaries**"""
    #return AbstractFilter(y, spec)
end

function morlet{T<:FFTW.fftwReal}(::FourierDomain{2},
        θ::T, aspectratio::T, center::T, den::T, N::Int, nPeriods::Int)
    halfN = N >> 1
    pstart = - ((nPeriods-1)>>1)
    pstop = (nPeriods-1)>>1 + iseven(nPeriods)
    ωstart = - halfN + pstart * N
    ωstop = halfN + pstop * N - 1
    ωrange = ωstart:ωstop
    nΩs = length(ωrange)
    ωrows = repmat(-ωrange, 1, nΩs)
    ωcols = repmat(ωrange', nΩs, 1)
    cosθ = cos(θ)
    sinθ = sin(θ)
    ωrows_θ = ωrows * cosθ + ωcols * sinθ
    ωcols_θ = -ωrows * sinθ + ωrows * cosθ
end
