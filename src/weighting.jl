abstract AbstractWeighting

immutable EqualWeighting <: AbstractWeighting end
weight_frequencies(::EqualWeighting, ξs) = ones(ξs)

immutable LoudnessWeighting <: AbstractWeighting
    samplerate::Int
end

function weight_frequencies(weighting::LoudnessWeighting, ξs)
    freqs = weighting.samplerate * ξs
    freqs2 = freqs .* freqs
    numerator = 12200.0 * freqs2
    sqden1 = (freqs2 + 20.6 * 20.6)
    sqden2 = sqrt.((freqs2 + 107.7 * 107.7) .* (freqs2 + 737.9 * 737.9))
    sqden3 = (freqs2 + 12200.0 * 12200.0)
    sqdenominator = sqden1 .* sqden2 .* sqden3
    denominator = sqrt.(sqdenominator)
    return numerator ./ denominator
end
