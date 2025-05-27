function get_sampled_tree(event_log::Vector{<:AbstractEvent})::Tree
    n_leaves = n_sampled(event_log)
    nodes = Vector{AbstractNode}(undef, 2*n_leaves)
    branch_lengths = Vector{Float64}(undef, 2*n_leaves-1)
    node_map = Dict{Int, Int}()
    node_id = 0
    for event in reverse(event_log)
        if event isa Sampling
            node_id += 1
            nodes[node_id] = Node(node_id, nothing, nothing, event.time)
            node_map[event.host] = node_id
        elseif event isa Transmission && haskey(node_map, event.infectee)
            if haskey(node_map, event.infector)
                node_id += 1
                left, right = node_map[event.infector], node_map[event.infectee]
                nodes[node_id] = Node(node_id, left, right, event.time)
                branch_lengths[left] = nodes[left].time - event.time
                branch_lengths[right] = nodes[right].time - event.time
                node_map[event.infector] = node_id
            else
                node_map[event.infector] = node_map[event.infectee]
            end
        elseif event isa Seed && haskey(node_map, event.host)
            node_id += 1
            left = node_map[event.host]
            nodes[node_id] = Node(node_id, left, nothing, event.time)
            branch_lengths[left] = nodes[left].time - event.time
        end
    end
    return Tree(nodes, branch_lengths)
end


@forward Simulation.event_log get_sampled_tree


function get_tree_stats(tree::Tree)
    root_height = tree[end].time
    root_to_tip_distances = [node.time - root_height for node in tree if isleaf(node)]
    n_children = Dict{Int, Int}()
    ladder_length = Dict{Int, Int}()
    n_leaves = 0
    n_cherries = 0
    n_ladder_nodes = 0

    colless_index = 0
    max_ladder_length = 0
    for node in tree
        if isleaf(node)
            n_leaves += 1
            n_children[node.id] = 1
            ladder_length[node.id] = 0
        elseif isbinary(node)
            n_children[node.id] = n_children[node.left] + n_children[node.right]
            colless_index += abs(n_children[node.left] - n_children[node.right])
            if isleaf(tree[node.left]) && isleaf(tree[node.right])
                n_cherries += 1
                ladder_length[node.id] = 0
            elseif isleaf(tree[node.left])
                n_ladder_nodes += 1
                ladder_length[node.id] = ladder_length[node.right] + 1
            elseif isleaf(tree[node.right])
                n_ladder_nodes += 1
                ladder_length[node.id] = ladder_length[node.left] + 1
            else
                ladder_length[node.id] = 0
            end
            if ladder_length[node.id] > max_ladder_length
                max_ladder_length = ladder_length[node.id]
            end
        end
    end

    tree_height = tree[1].time - tree[end].time
    # colless_index /= (n_leaves - 1) * (n_leaves - 2) / 2
    sackin_index = sum(root_to_tip_distances) / n_leaves
    return tree_height, colless_index, sackin_index, n_leaves, n_cherries, n_ladder_nodes, max_ladder_length
end