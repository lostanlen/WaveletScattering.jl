abstract AbstractPointGroup
immutable TrivialGroup <: AbstractPointGroup end
immutable ReflectionGroup <: AbstractPointGroup end
immutable RotationGroup <: AbstractPointGroup
    nOrientations::Int
end

typealias LineGroups Union{TrivialGroup,ReflectionGroup}

get_nOrientations(group::TrivialGroup) = 1
get_nOrientations(group::ReflectionGroup) = 2
get_nOrientations(group::RotationGroup) = group.nOrientations
