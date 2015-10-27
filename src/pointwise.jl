abstract AbstractPointwise
immutable Modulus <: AbstractPointwise end

map!(ρ::Modulus, blob_in::AbstractNode, blob_out::AbstractNode) =
    map!(abs, blob_in.data, blob_out.data)
