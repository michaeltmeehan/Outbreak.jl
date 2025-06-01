function get_sampled_tree(event_log::Vector{<:AbstractEvent})::Tree
    n_leaves = n_sampled(event_log)
    nodes = Vector{AbstractNode}(undef, 2*n_leaves)
    n_leaves == 0 && return Tree(nodes, Float64[])
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
    
    depth, width = get_tree_size(tree)

    tree_height = tree[1].time - tree[end].time
    colless_index /= (n_leaves - 1) * (n_leaves - 2) / 2
    sackin_index = sum(root_to_tip_distances) / n_leaves
    ladder_length = max_ladder_length / n_leaves
    il_portion = n_ladder_nodes / n_leaves
    max_width_on_depth = maximum(width) / maximum(depth)
    max_width_difference = maximum(diff(width))
    return (; tree_height = tree_height,
            colless_index = colless_index,
            sackin_index = sackin_index,
            n_leaves = n_leaves,
            n_cherries = n_cherries,
            n_ladder_nodes = n_ladder_nodes,
            max_ladder_length = max_ladder_length,
            ladder_length = ladder_length,
            il_portion = il_portion,
            max_width_on_depth = max_width_on_depth,
            max_width_difference = max_width_difference)
end


function get_tree_size(tree::Tree)
    depth = Vector{Int}(undef, length(tree))
    for node in reverse(tree)
        if isroot(node)
            depth[node.id] = 0
            depth[node.left] = 1
        elseif isbinary(node)
            depth[node.left] = depth[node.right] = depth[node.id] + 1
        end
    end
    width = _tabulate(depth)
    return depth, width
end


function _tabulate(v::Vector{Int})::Vector{Int}
    out = fill(0, maximum(v) + 1)
    for i in v
        out[i+1] += 1
    end
    return out
end


function _vectorize(d::Dict{Int, Int})::Vector{Int}
    out = fill(0, maximum(keys(d)) + 1)
    for (k, v) in d
        out[k+1] = v
    end
    return out
end


function get_ltt(tree::Tree)::Tuple{Vector{Float64}, Vector{Int}}
    t = Vector{Float64}(undef, length(tree))
    n = fill(0, length(tree))
    for (i, node) in enumerate(reverse(tree))
        t[i] = node.time
        if isleaf(node)
            n[i] = n[i-1] - 1
        elseif isbinary(node)
            n[i] = n[i-1] + 1
        elseif isroot(node)
            n[i] = i > 1 ? n[i-1] + 1 : 1
        end
    end
    return t, n
end