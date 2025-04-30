
function filter_event_log(event_log::Vector{<: EpiSim.EpiEvent.AbstractEpiEvent})
    filtered_log = Vector{EpiSim.EpiEvent.AbstractEpiEvent}()
    ancestors = Set{Int}()
    for event in reverse(event_log)
        if event isa EpiSim.Sampling
            push!(filtered_log, event)
            push!(ancestors, event.host)
        elseif event isa EpiSim.Transmission && event.infectee in ancestors
                push!(filtered_log, event)
                push!(ancestors, event.infector)
        elseif event isa EpiSim.Seed && event.host in ancestors
                push!(filtered_log, event)
                push!(ancestors, event.host)
        end
    end
    return reverse(filtered_log)
end


function add_leaf!(tree::Vector{<:AbstractNode}, event::EpiSim.Sampling, node_map::Dict{Int, Int}, node_id::Int)
    node_id += 1
    push!(tree, LeafNode(node_id, event.time))
    node_map[event.host] = node_id
    return node_id
end


function add_binary!(tree::Vector{<:AbstractNode}, event::EpiSim.Transmission, node_map::Dict{Int, Int}, node_id::Int)
    if haskey(node_map, event.infector)
        node_id += 1
        push!(tree, BinaryNode(node_id, node_map[event.infector], node_map[event.infectee], event.time))
        node_map[event.infector] = node_id
    else
        node_map[event.infector] = node_map[event.infectee]
    end
    return node_id
end


function add_root!(tree::Vector{<:AbstractNode}, event::EpiSim.Seed, node_map::Dict{Int, Int}, node_id::Int)
    node_id += 1
    push!(tree, RootNode(node_id, nothing, node_map[event.host], event.time))
    return node_id
end


function get_sampled_tree(event_log::Vector{<:EpiSim.EpiEvent.AbstractEpiEvent})::Vector{<:AbstractNode}
    tree = Vector{AbstractNode}()
    node_map = Dict{Int, Int}()
    node_id = 0
    for event in reverse(event_log)
        if event isa EpiSim.Sampling
            node_id = add_leaf!(tree, event, node_map, node_id)
        elseif event isa EpiSim.Transmission && haskey(node_map, event.infectee)
            node_id = add_binary!(tree, event, node_map, node_id)
        elseif event isa EpiSim.Seed && haskey(node_map, event.host)
            node_id = add_root!(tree, event, node_map, node_id)
        end
    end
    return tree
end