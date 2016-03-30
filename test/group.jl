using Base.Test
import WaveletScattering: get_nOrientations, ReflectionGroup,
    RotationGroup, TrivialGroup

@test get_nOrientations(TrivialGroup()) == 1
@test get_nOrientations(ReflectionGroup()) == 2
@test get_nOrientations(RotationGroup(4)) == 4
@test get_nOrientations(RotationGroup(8)) == 8
