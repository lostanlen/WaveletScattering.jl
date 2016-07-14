for T in subtypes(AbstractPointwise)
    @eval begin
        function (ρ::$T)(node::AbstractNode)
            return Node(ρ(node.data), node.ranges)
        end
    end
end
