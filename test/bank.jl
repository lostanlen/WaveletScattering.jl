import WaveletScattering: kthrange, PathKey, Literal

r = kthrange([:time], rand(2), 1)
@test r == (PathKey(Literal[Literal(:time,1)])=>1:1:2)
