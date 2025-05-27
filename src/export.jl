format_subtree(node_id::Int, branch_length::Float64, node_strings::Dict{Int, String}) =
    "$(node_strings[node_id]):$(round(branch_length; sigdigits=4))"


function build_newick(tree::Tree{<:AbstractNode})::String
    node_strings = Dict{Int, String}()

    for node in tree
        if isleaf(node)
            node_strings[node.id] = string(node.id)

        elseif isbinary(node)
            left_subtree = format_subtree(node.left, tree.branch_lengths[node.left], node_strings)
            right_subtree = format_subtree(node.right, tree.branch_lengths[node.right], node_strings)
            node_strings[node.id] = "($left_subtree,$right_subtree)"

        elseif isroot(node)
            parts = String[]
            !isnothing(node.left) && push!(parts, format_subtree(node.left, tree.branch_lengths[node.left], node_strings))
            !isnothing(node.right) && push!(parts, format_subtree(node.right, tree.branch_lengths[node.right], node_strings))
            node_strings[node.id] = "(" * join(parts, ",") * ")"
        end
    end

    # Support forests (multiple roots)
    roots = [node.id for node in tree if isroot(node)]
    top_strings = [node_strings[id] for id in roots]
    return join(top_strings, "*") * ";"
end


function write_nexus(filename::String,
                     tree::Tree{<:AbstractNode},
                     alignment::Vector{Sequence})

    open(filename, "w") do io
        # Write alignment block using SeqSim
        SeqSim.write_nexus(io, alignment)

        # Write trees block
        println(io, "\nBEGIN TREES;")
        println(io, "    TREE tree_1 = [&R] ", build_newick(tree))
        println(io, "END;")
    end
end
