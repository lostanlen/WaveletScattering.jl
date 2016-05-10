using Base.Test
# blob.jl
import WaveletScattering: ScatteredBlob
# node.jl
import WaveletScattering: Node
# path.jl
import WaveletScattering: Path, PathKey

# show
io = IOBuffer()
node = Node([1.0, 2.0], (PathKey(:time)=>1:1:2,))
blob = ScatteredBlob(DataStructures.SortedDict((Path() => node,)))
show(io, blob)
@test takebuf_string(io) == "ScatteredBlob(1 node)"
blob = ScatteredBlob(DataStructures.SortedDict(
    (Path((:γ, :time) => 0) => node, Path((:γ, :time) => 1) => node)))
show(io, blob)
@test takebuf_string(io) == "ScatteredBlob(2 nodes)"
