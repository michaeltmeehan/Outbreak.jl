format_subtree(child_id::Int, parent::AbstractNode, tree::Vector{<:AbstractNode}, node_strings::Dict{Int, String}) =
    "$(node_strings[child_id]):$(round(get_branch_length(parent, child_id, tree); sigdigits=4))"


function build_newick(tree::Vector{<:AbstractNode})::String
    @assert issorted(tree) "Tree is not sorted in descending order of time."
    node_strings = Dict{Int, String}()

    for node in tree
        if isleaf(node)
            node_strings[node.id] = string(node.id)

        elseif isbinary(node)
            left_subtree = format_subtree(node.left, node, tree, node_strings)
            right_subtree = format_subtree(node.right, node, tree, node_strings)
            node_strings[node.id] = "($left_subtree,$right_subtree)"

        elseif isroot(node)
            parts = String[]
            if node.left !== nothing
                push!(parts, format_subtree(node.left, node, tree, node_strings))
            end
            push!(parts, format_subtree(node.right, node, tree, node_strings))
            node_strings[node.id] = "(" * join(parts, ",") * ")"
        end
    end

    # Support forests (multiple roots)
    roots = get_root_id(tree)
    top_strings = [node_strings[id] for id in roots]
    return join(top_strings, "*") * ";"
end


function write_nexus(filename::String,
                     tree::Vector{<:AbstractNode},
                     alignment::Vector{Sequence})
    @assert issorted(tree) "Tree must be sorted in descending order of time."

    open(filename, "w") do io
        # Write alignment block using SeqSim
        SeqSim.write_nexus(io, alignment)

        # Write trees block
        println(io, "\nBEGIN TREES;")
        println(io, "    TREE tree_1 = [&R] ", build_newick(tree))
        println(io, "END;")
    end
end
