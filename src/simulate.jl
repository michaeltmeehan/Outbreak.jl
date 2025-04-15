function simulate_alignment(rng::AbstractRNG,
                            event_log::Vector{<:EpiSim.AbstractEpiEvent},
                            site_model::SeqSim.SiteModel;
                            root_seq::Union{Nothing, SeqSim.Sequence}=nothing)::Vector{SeqSim.Sequence}

    # Initialize an empty sequence tree
    sequence_tree = Dict{Int, SeqSim.Sequence}()

    # Initialize the root sequence in the tree
    alignment = Vector{SeqSim.Sequence}()

    # Unpack site_model
    sequence_length = site_model.sequence_length
    μ = site_model.μ
    variable_sites = site_model.variable_sites
    gamma_category_count = site_model.gamma_category_count
    λ = site_model.substitution_model.λ
    V = site_model.substitution_model.V
    V⁻¹ = site_model.substitution_model.V⁻¹
    
    # Initialize transition weights
    transition_weights = [zeros(4,4) for _ in 1:gamma_category_count]

    for event in event_log
        if event isa EpiSim.Seed
            sequence_tree[event.host] = isnothing(root_seq) ? SeqSim.rand_seq(rng, sequence_length, taxon=event.host, time=event.time) : root_seq

        elseif event isa EpiSim.Transmission
            # Get the parent sequence from the tree
            parent_sequence = sequence_tree[event.infector]

            # calculate the time since update
            Δt = event.time - parent_sequence.time

            # Update parent sequence and pass it to child
            SeqSim.update_sequence!(rng, parent_sequence.value, transition_weights, Δt, μ, variable_sites, λ, V, V⁻¹)
            sequence_tree[event.infectee] = SeqSim.Sequence(event.infectee, copy(parent_sequence.value), event.time)

        elseif event isa EpiSim.Sampling
            # Get the sampled sequence from the tree
            sampled_sequence = sequence_tree[event.host]

            # Calculate the time since last update
            Δt = event.time - sampled_sequence.time

            # Update sequence
            SeqSim.update_sequence!(rng, sampled_sequence.value, transition_weights, Δt, μ, variable_sites, λ, V, V⁻¹)

            # Add the sampled sequence to the alignment
            push!(alignment, sampled_sequence)
        end
    end
    return alignment
end