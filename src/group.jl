abstract AbstractPointGroup
immutable TrivialGroup <: AbstractPointGroup end
immutable ReflectionGroup <: AbstractPointGroup end
immutable RotationGroup <: AbstractPointGroup end

typealias LineGroups Union{TrivialGroup,ReflectionGroup}
