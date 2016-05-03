for T in subtypes(AbstractPointwise)
    @eval begin
        function call{NODE<:AbstractNode}(ρ::$T, node::NODE)
            return Node(ρ(node.data), node.ranges)
        end
    end
end
