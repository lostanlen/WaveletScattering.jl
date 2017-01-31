using Base.Test
import WaveletScattering: get_n_orientations, ReflectionGroup,
    RotationGroup, TrivialGroup

@test get_n_orientations(TrivialGroup()) == 1
@test get_n_orientations(ReflectionGroup()) == 2
@test get_n_orientations(RotationGroup(4)) == 4
@test get_n_orientations(RotationGroup(8)) == 8
