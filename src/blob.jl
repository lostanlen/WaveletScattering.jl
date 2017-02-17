# ScatteredBlob
immutable ScatteredBlob{NODE<:AbstractNode} <: Mocha.Blob
    nodes::DataStructures.SortedDict{Path,NODE}
end

function Base.show(io::IO, blob::ScatteredBlob)
    n_nodes = length(blob.nodes)
    plural = n_nodes > 1 ? "s" : ""
    print(io, "ScatteredBlob(", n_nodes, " node", plural, ")")
end

function forward!(
        backend,
        blob_out::ScatteredBlob,
        bank::Bank1D,
        blob_in::ScatteredBlob)
    pathdepth(blob_in, bank.behavior.pathkey)
    map(node -> pathdepth(bank.behavior.pathkey, keys(blob_in.nodes)))
    γkey = cons(Literal(:γ, 1), bank.behavior.pathkey)
    for j in bank.behavior.j_range
        for χ in 0:(bank.spec.n_filters_per_octave-1)
            ψ = bank.ψs[1 + θ, 1 + χ, 1 + j]
            for (path_in, node_in) in input.nodes
                path_out = copy(path_in)
                path_out[γkey] = γ
                transform!(blob[path_out], blob_in[path_in], ψ)
            end
        end
    end
end
