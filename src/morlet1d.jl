"""Computes `gauss{T<:Real}(ω, den::T) = exp(- ω*ω/den)`.
The Gaussian bell curve is defined as `gauss(ω) = exp(- ω² / 2σ²)`.
For performance reasons, we memoize the denominator `2σ²`,
which is computed only once in the caller `morlet`.
Also note that the exponentiation `ω²` is replaced by the explicit
product `ω*ω`.
"""
gauss{T<:Real}(ω, den::T) = @fastmath Base.convert(T, exp(- ω*ω/den))::T

"""Computes a one-dimensional Morlet wavelet in the Fourier domain.
A Morlet wavelet of center frequency `ξ` and of variance `σ` looks almost like
a Gaussian bell curve. To ensure that the wavelet has a vanishing moment, we
substract a corrective term around the zeroth frequency."""
function AbstractFilter{T<:FFTW.fftwReal,G<:LineGroups}(ψmeta::ΨMeta1D,
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
    """2. **Number of periods**"""
    halfsupport = sqrt(den * log(inv(spec.ɛ)))
    firstω = max(center - halfsupport, -5N/2)
    lastω = min(center + halfsupport, 5N/2 - 1)
    nPeriods = 1 + ceil(Int, (lastω-halfN) / N)
    """3. **Call to morlet**"""
    y = morlet(FourierDomain(1), center, den, N, nPeriods)
    """4. **Trimming to true support boundaries**"""
    return AbstractFilter(y, spec)
end

function AbstractFilter{T<:FFTW.fftwReal,G<:LineGroups}(ϕmeta::ΦMeta,
        spec::Spec1D{T,FourierDomain{1},G,Morlet})
    """1. **Gaussian denominator `den = 2σ²`**"""
    @inbounds log2_length = spec.log2_size
    N = 1 << log2_length
    halfN = N >> 1
    bw = N * T(ϕmeta.bandwidth)
    den = @fastmath bw * bw / T(2.0 * log(2.0))
    """2. **Call to morlet**"""
    lastω = round(Int, min(bw * sqrt(2.0 / spec.ɛ), halfN - 1))
    leg = T[ gauss(ω, den) for ω in 1:lastω ]
    """3. **Trimming to true support boundaries**"""
    leg = leg[1:findlast(leg .> spec.ɛ)]
    return FourierSymmetric1DFilter(leg, one(T))
end

function morlet{T<:FFTW.fftwReal}(::FourierDomain{1},
        center::T, den::T, N::Int, nPeriods::Int)
    halfN = N >> 1
    pstart = - ((nPeriods-1)>>1)
    pstop = (nPeriods-1)>>1 + iseven(nPeriods)
    ωstart = - halfN + pstart * N
    ωstop = halfN + pstop * N - 1
    @inbounds begin
        gauss_center = T[ gauss(ω-center, den) for ω in ωstart:ωstop]
        gauss_0 = T[ gauss(ω, den)
            for ω in (ωstart + pstart*N):(ωstop + pstop*N) ]
        corrective_gaussians = T[ gauss_0[1 + ω + p*N]
            for ω in 0:(N*nPeriods-1), p in 0:(nPeriods-1) ]
    end
    b = T[ gauss(p*N - center, den) for p in pstart:pstop ]
    A = T[ gauss((q-p)*N, den)
        for p in 0:(nPeriods-1), q in 0:(nPeriods-1) ]
    corrective_factors = A \ b
    y = gauss_center - corrective_gaussians * corrective_factors
    y = reshape(y, N, nPeriods)
    y = squeeze(sum(y, 2), 2)
end
