"""Three variables come into play in a wavelet filterbank:
* the spatial variable `V`, which may be multi-dimensional ;
* the scale variable `γ::V`, which is one-dimensional and non-negative ;
* the isometric variable `θ::V`, which belongs to a ""point group"".
A point group is an algebraic group of transformations which have a single
fixed point. Reflections of the real line and rotations of the plane are
examples of such point groups.
Conversely, translations and scalings are not point groups.
"""
abstract AbstractPointGroup

"""The `TrivialGroup` is a finite group of a single neutral element."""
immutable TrivialGroup <: AbstractPointGroup end

"""The `ReflectionGroup` is a finite group of two elements, which represents
the reflections of the real line."""
immutable ReflectionGroup <: AbstractPointGroup end

"""The `RotationGroup` is a continuous group, which represents the rotations of
the plane. Its orientations are discretized uniformly."""
immutable RotationGroup <: AbstractPointGroup
    n_orientations::Int
end

typealias LineGroups Union{TrivialGroup,ReflectionGroup}
typealias PlaneGroups Union{TrivialGroup,RotationGroup}

"""Returns the number of orientations associated to a point group.
Equal to `1` for `TrivialGroup`, `2` for `ReflectionGroup`,
and `group.n_orientations` for `RotationGroup`."""
get_n_orientations(group::TrivialGroup) = 1
get_n_orientations(group::ReflectionGroup) = 2
get_n_orientations(group::RotationGroup) = group.n_orientations
