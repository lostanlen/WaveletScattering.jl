using Base.Test
# behavior.jl
import WaveletScattering: Behavior
# bank.jl
import WaveletScattering: renormalize!
# filter.jl
import WaveletScattering: AbstractFilter
# path.jl
import WaveletScattering: PathKey
# spec1d.jl
import WaveletScattering: Spec1D
# weighting.jl
import WaveletScattering: EqualWeighting

pathkey = PathKey(:time)
spec = Spec1D()
(nΘs, nΧs, nJs) = size(spec.ψmetas)
domaintype = typeof(spec.domain)
ψs = Array(AbstractFilter{spec.signaltype,domaintype}, (nΘs, nΧs, nJs))
ψs[1, :, :] =
    pmap(AbstractFilter, spec.ψmetas[1, :, :], fill(spec, nΧs, nJs))
(nΘs > 1) && spin!(ψs)
ϕ = AbstractFilter(spec.ϕmeta, spec)
renormalize!(ϕ, ψs, spec)

is_ϕ_applied = false
j_range = 0:(spec.n_octaves-1)
weighting = EqualWeighting()

for log2_oversampling = 0:5
    for max_log2_stride = 0:(spec.n_octaves-1)
        behavior = Behavior(ϕ, ψs, spec, is_ϕ_applied, j_range,
            log2_oversampling, max_log2_stride, pathkey, weighting)
        @test all(behavior.ψ_log2_samplings .<= 0)
        @test all(behavior.ψ_log2_samplings .>= -max_log2_stride)
    end
end

behavior = Behavior(ϕ, ψs, spec, is_ϕ_applied, j_range,
    0, (spec.n_octaves-1), pathkey, weighting)
@test all(behavior.weights .== 1.0)
