
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


struct Node
    left::Union{Nothing, Int}
    right::Union{Nothing, Int}
    time::Float64
end


function get_sampled_tree(event_log::Vector{<:EpiSim.EpiEvent.AbstractEpiEvent})::Vector{Node}
    ancestors = Set{Int}()
    tree = Vector{Node}()
    node_map = Dict{Int, Int}()
    for event in reverse(event_log)
        if event isa EpiSim.Sampling
            push!(tree, Node(nothing, nothing, event.time))
            node_map[event.host] = length(tree)
            push!(ancestors, event.host)
        elseif event isa EpiSim.Transmission && event.infectee ∈ ancestors
            pop!(ancestors, event.infectee)
            if event.infector in ancestors
                push!(tree, Node(node_map[event.infector], node_map[event.infectee], event.time))
                node_map[event.infector] = length(tree)
            else
                node_map[event.infector] = node_map[event.infectee]
                push!(ancestors, event.infector)
            end
        elseif event isa EpiSim.Seed && event.host ∈ ancestors
            push!(tree, Node(nothing, node_map[event.host], event.time))
            push!(ancestors, event.host)
        end
    end
    return tree
end