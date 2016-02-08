abstract AbstractWeighting

immutable EqualWeighting <: AbstractWeighting
weight_frequencies(::AbstractWeighting, ξs) = ones(ξs)

function weight_frequencies(weighting::EqualWeighting)

immutable LoudnessWeighting
    samplerate::Int
end

function weight_frequencies(weighting::LoudnessWeighting, ξs)
    freqs = weighting.samplerate * ξs
    freqs2 = freqs .* freqs
    freqs4 = freqs2 .* freqs2
    numerator = 12200.0 * 12200.0 * freqs
    den1 = (freqs2 + 20.6 * 20.6)
    den2 = sqrt((freqs2 + 107.7 * 107.7) * (freqs2 + 737.9 * 737.9))
    den3 = (freqs2 + 12200.0 * 12200.0)
    denominator = den1 .* den2 .* den3
    return convert(numerator ./ denominator
end
