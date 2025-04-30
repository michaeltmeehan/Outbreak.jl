function build_newick(tree::Vector{Node})::String
    node_strings = Dict{Int, String}()
    
    for node in reverse(tree)  # start from root (first) down to tips
        if isleaf(node)
            label = isnothing(node.taxon) ? "taxon$(node.id)" : string(node.taxon)
            node_strings[node.id] = label
        else
            left_str = node_strings[node.left]
            right_str = node_strings[node.right]
            node_strings[node.id] = "($left_str,$right_str)"
        end
    end
    
    return node_strings[tree[1].id] * ";"  # tree[1] is root due to sorting
end
