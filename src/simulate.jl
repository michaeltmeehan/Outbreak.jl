function propagate_to_child!(rng::AbstractRNG, tree::Tree{<:AbstractNode}, parent_node::AbstractNode, child_id::Union{Nothing, Int}, sequence::Vector{UInt8},
                          sequence_tree::Dict{Int, Vector{UInt8}}, prop::SeqSim.SequencePropagator)
    Δt = get_branch_length(parent_node, child_id, tree)
    @assert Δt ≥ 0 "Negative branch length encountered between node $(parent_node.id) and $child_id."
    sequence_tree[child_id] = prop(rng, sequence, Δt)
end


function propagate_to_children!(rng::AbstractRNG, tree::Tree{<:AbstractNode}, parent_node::AbstractNode, sequence::Vector{UInt8}, sequence_tree::Dict{Int, Vector{UInt8}}, prop::SequencePropagator)
    !isnothing(parent_node.left) && propagate_to_child!(rng, tree, parent_node, parent_node.left, sequence, sequence_tree, prop)
    !isnothing(parent_node.right) && propagate_to_child!(rng, tree, parent_node, parent_node.right, sequence, sequence_tree, prop)
end


function simulate_alignment(rng::AbstractRNG,
                            tree::Tree{<:AbstractNode},
                            prop::SequencePropagator,
                            root_seq::Union{Nothing, Vector{UInt8}})::Vector{Sequence}
    sequence_tree = Dict{Int, Vector{UInt8}}()
    alignment = Vector{Sequence}()

    # Main propagation loop
    for node in reverse(tree)
        sequence = get!(sequence_tree, node.id) do
            isnothing(root_seq) ? SeqSim.rand_seq_int(rng, prop.site_model) : root_seq
        end
        if isleaf(node)
            push!(alignment, Sequence(sequence_tree[node.id], taxon=node.id, time=node.time))
        else
            propagate_to_children!(rng, tree, node, sequence, sequence_tree, prop)
        end
    end
    return alignment
end


function simulate_alignment(tree::Tree{<:AbstractNode},
                            prop::SequencePropagator;
                            root_seq::Union{Nothing, Vector{UInt8}}=nothing)::Alignment
    return simulate_alignment(Random.GLOBAL_RNG, tree, prop, root_seq)
end