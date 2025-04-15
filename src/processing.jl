
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