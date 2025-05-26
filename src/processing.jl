
# function filter_event_log(event_log::Vector{<: EpiSim.EpiEvent.AbstractEpiEvent})
#     filtered_log = Vector{EpiSim.EpiEvent.AbstractEpiEvent}()
#     ancestors = Set{Int}()
#     for event in reverse(event_log)
#         if event isa EpiSim.Sampling
#             push!(filtered_log, event)
#             push!(ancestors, event.host)
#         elseif event isa EpiSim.Transmission && event.infectee in ancestors
#                 push!(filtered_log, event)
#                 push!(ancestors, event.infector)
#         elseif event isa EpiSim.Seed && event.host in ancestors
#                 push!(filtered_log, event)
#                 push!(ancestors, event.host)
#         end
#     end
#     return reverse(filtered_log)
# end


# function add_leaf!(tree::Vector{<:AbstractNode}, event::Sampling, node_map::Dict{Int, Int}, node_id::Int)
#     node_id += 1
#     push!(tree, LeafNode(node_id, event.time))
#     node_map[event.host] = node_id
#     return node_id
# end


# # TODO: Perhaps rename this because a binary node is not always added (i.e., the else branch)
# function add_binary!(tree::Vector{<:AbstractNode}, event::Transmission, node_map::Dict{Int, Int}, node_id::Int)
#     if haskey(node_map, event.infector)
#         node_id += 1
#         push!(tree, BinaryNode(node_id, node_map[event.infector], node_map[event.infectee], event.time))
#         node_map[event.infector] = node_id
#     else
#         node_map[event.infector] = node_map[event.infectee]
#     end
#     return node_id
# end


# function add_root!(tree::Vector{<:AbstractNode}, event::Seed, node_map::Dict{Int, Int}, node_id::Int)
#     node_id += 1
#     push!(tree, RootNode(node_id, nothing, node_map[event.host], event.time))
#     return node_id
# end


# function get_sampled_tree(event_log::Vector{<:AbstractEvent})::Vector{<:AbstractNode}
#     tree = Vector{AbstractNode}()
#     node_map = Dict{Int, Int}()
#     node_id = 0
#     for event in reverse(event_log)
#         if event isa EpiSim.Sampling
#             node_id = add_leaf!(tree, event, node_map, node_id)
#         elseif event isa EpiSim.Transmission && haskey(node_map, event.infectee)
#             node_id = add_binary!(tree, event, node_map, node_id)
#         elseif event isa EpiSim.Seed && haskey(node_map, event.host)
#             node_id = add_root!(tree, event, node_map, node_id)
#         end
#     end
#     return tree
# end


function get_sampled_tree(event_log::Vector{<:AbstractEvent})::Tree
    tree = Vector{AbstractNode}()
    node_map = Dict{Int, Int}()
    for event in reverse(event_log)
        if event isa Sampling
            node_id = length(tree) + 1
            push!(tree, Node(node_id, nothing, nothing, event.time))
            node_map[event.host] = node_id
        elseif event isa Transmission && haskey(node_map, event.infectee)
            if haskey(node_map, event.infector)
                node_id = length(tree) + 1
                push!(tree, Node(node_id, node_map[event.infector], node_map[event.infectee], event.time))
                node_map[event.infector] = node_id
            else
                node_map[event.infector] = node_map[event.infectee]
            end
        elseif event isa Seed && haskey(node_map, event.host)
            node_id = length(tree) + 1
            push!(tree, Node(node_id, nothing, node_map[event.host], event.time))
        end
    end
    return Tree(tree)
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