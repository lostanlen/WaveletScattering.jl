import WaveletScattering: AbstractDomain, FourierDomain, LineDomains,
    PlaneDomains, SpatialDomain

# FourierDomain
@test issubtype(FourierDomain, AbstractDomain)
@test FourierDomain(1).dim == Val{1}()
@test FourierDomain(2).dim == Val{2}()
@test_throws AssertionError FourierDomain(0)
@test_throws AssertionError FourierDomain(-1)
@test_throws MethodError FourierDomain(0.5)

# SpatialDomain
@test issubtype(SpatialDomain, AbstractDomain)
@test SpatialDomain(1).dim == Val{1}()
@test SpatialDomain(2).dim == Val{2}()
@test_throws AssertionError SpatialDomain(0)
@test_throws AssertionError SpatialDomain(-1)
@test_throws MethodError FourierDomain(0.5)
