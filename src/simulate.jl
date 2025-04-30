function propagate_to_child!(rng::AbstractRNG, tree::Vector{Node}, parent_node::Node, child_id::Union{Nothing, Int}, sequence::Vector{UInt8},
                          sequence_tree::Dict{Int, Vector{UInt8}}, prop::SeqSim.SequencePropagator)
    isnothing(child_id) && return
    Δt = tree[child_id].time - parent_node.time
    @assert Δt ≥ 0 "Negative branch length encountered between node $(parent_node.id) and $child_id."
    sequence_tree[child_id] = prop(rng, sequence, Δt)
end


function simulate_alignment(rng::AbstractRNG,
                            tree::Vector{Node},
                            prop::SeqSim.SequencePropagator,
                            root_seq::Union{Nothing, Vector{UInt8}})::Vector{SeqSim.Sequence}
    @assert issorted(tree) "Tree is not sorted in descending order of time."
    sequence_tree = Dict{Int, Vector{UInt8}}()
    alignment = Vector{SeqSim.Sequence}()

    # Main simulation loop
    for node in reverse(tree)
        sequence = get!(sequence_tree, node.id) do
            isnothing(root_seq) ? SeqSim.rand_seq_int(rng, prop.site_model) : root_seq
        end

        propagate_to_child!(rng, tree, node, node.left, sequence, sequence_tree, prop)
        propagate_to_child!(rng, tree, node, node.right, sequence, sequence_tree, prop)
        
        isleaf(node) && push!(alignment, Sequence(sequence_tree[node.id], taxon=node.id, time=node.time))
    end
    return alignment
end


function simulate_alignment(tree::Vector{Node},
                            prop::SeqSim.SequencePropagator;
                            root_seq::Union{Nothing, Vector{UInt8}}=nothing)::Vector{SeqSim.Sequence}
    return simulate_alignment(Random.GLOBAL_RNG, tree, prop, root_seq)
end