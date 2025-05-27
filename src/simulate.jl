function simulate_alignment(rng::AbstractRNG,
                            tree::Tree{<:AbstractNode},
                            seqprop::SequencePropagator,
                            root_seq::Union{Nothing, Vector{UInt8}})::Vector{Sequence}
    sequences = Dict{Int, Vector{UInt8}}()
    alignment = Vector{Sequence}()

    # Main propagation loop
    for node in reverse(tree)
        sequence = get!(sequences, node.id) do
            isnothing(root_seq) ? SeqSim.rand_seq_int(rng, seqprop.site_model) : root_seq
        end
        if isleaf(node)
            push!(alignment, Sequence(sequences[node.id], taxon=node.id, time=node.time))
        elseif isbinary(node)
            sequences[node.left] = seqprop(rng, sequence, tree.branch_lengths[node.left])
            sequences[node.right] = seqprop(rng, sequence, tree.branch_lengths[node.right])
        end
    end
    return alignment
end


function simulate_alignment(tree::Tree{<:AbstractNode},
                            seqprop::SequencePropagator;
                            root_seq::Union{Nothing, Vector{UInt8}}=nothing)::Alignment
    return simulate_alignment(Random.GLOBAL_RNG, tree, seqprop, root_seq)
end