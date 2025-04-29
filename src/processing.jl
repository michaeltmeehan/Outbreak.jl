
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


function get_sampled_tree(event_log::Vector{<:EpiSim.EpiEvent.AbstractEpiEvent})::Vector{Node}
    tree = Vector{Node}()
    node_map = Dict{Int, Int}()
    num_nodes = 0
    for event in reverse(event_log)
        if event isa EpiSim.Sampling
            num_nodes += 1
            push!(tree, Node(num_nodes, nothing, nothing, event.time))
            node_map[event.host] = num_nodes
        elseif event isa EpiSim.Transmission && haskey(node_map, event.infectee)
            if haskey(node_map, event.infector)
                num_nodes += 1
                push!(tree, Node(num_nodes, node_map[event.infector], node_map[event.infectee], event.time))
                node_map[event.infector] = num_nodes
            else
                node_map[event.infector] = node_map[event.infectee]
            end
        elseif event isa EpiSim.Seed && haskey(node_map, event.host)
            num_nodes += 1
            push!(tree, Node(num_nodes, nothing, node_map[event.host], event.time))
        end
    end
    return tree
end