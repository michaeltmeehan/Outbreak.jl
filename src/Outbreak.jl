module Outbreak

using Random
using SeqSim
using TreeSim

export Outbreak,
       has_log,
       has_tree,
       has_alignment,
       with_tree,
       with_alignment,
       outbreak_tree,
       outbreak_alignment,
       simulate_outbreak_tree,
       simulate_outbreak_alignment

"""
    Outbreak(log[, tree[, aln]])

Thin orchestration bundle for genomic outbreak workflows.

The fields hold the underlying package objects directly. Missing workflow
stages are represented by `nothing`, so `Outbreak(log)`, `Outbreak(log, tree)`,
and `Outbreak(log, tree, aln)` are all valid partial states. A bundle may also
hold a downstream object without an upstream one, such as `Outbreak(nothing,
tree, aln)`, when orchestrating from a user-supplied tree.
"""
struct Outbreak{L,T,A}
    log::L
    tree::T
    aln::A
end

Outbreak(log) = Outbreak(log, nothing, nothing)
Outbreak(log, tree) = Outbreak(log, tree, nothing)

has_log(outbreak::Outbreak) = outbreak.log !== nothing
has_tree(outbreak::Outbreak) = outbreak.tree !== nothing
has_alignment(outbreak::Outbreak) = outbreak.aln !== nothing

function Base.show(io::IO, outbreak::Outbreak)
    print(io, "Outbreak(")
    _print_stage(io, "log", outbreak.log)
    print(io, ", ")
    _print_stage(io, "tree", outbreak.tree)
    print(io, ", ")
    _print_stage(io, "aln", outbreak.aln)
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", outbreak::Outbreak)
    println(io, "Outbreak")
    println(io, "  log  : ", _stage_summary(outbreak.log))
    println(io, "  tree : ", _stage_summary(outbreak.tree))
    print(io, "  aln  : ", _stage_summary(outbreak.aln))
end

function _print_stage(io::IO, name::AbstractString, value)
    print(io, name, "=")
    value === nothing ? print(io, "empty") : print(io, _stage_type(value))
    return nothing
end

_stage_summary(::Nothing) = "empty"
_stage_summary(value) = string("populated (", _stage_type(value), ")")
_stage_type(value) = string(nameof(typeof(value)))
_stage_type(value::AbstractVector) = "Vector"

"""
    with_tree(outbreak; validate=true)
    with_tree(log; validate=true)

Populate the tree stage from an event log using `TreeSim.tree_from_eventlog`.
The event-log-to-tree semantics live in TreeSim and its package extensions.
"""
function with_tree(outbreak::Outbreak; validate::Bool=true)
    tree = TreeSim.tree_from_eventlog(outbreak.log; validate)
    return Outbreak(outbreak.log, tree, outbreak.aln)
end

with_tree(log; validate::Bool=true) = with_tree(Outbreak(log); validate)

"""
    outbreak_tree(log; validate=true)

Return only the sampled-ancestry tree extracted by TreeSim.
"""
outbreak_tree(log; validate::Bool=true) = TreeSim.tree_from_eventlog(log; validate)

"""
    simulate_outbreak_tree(log; validate=true)

Run the event-log to sampled-ancestry tree stage and return
`Outbreak(log, tree, nothing)`.
"""
function simulate_outbreak_tree(log; validate::Bool=true)
    return with_tree(log; validate)
end

"""
    with_alignment(rng, outbreak, site_model)
    with_alignment(outbreak, site_model)
    with_alignment(tree, site_model)

Populate the alignment stage using `SeqSim.simulate_alignment`. The tree-driven
sequence semantics live in SeqSim and its package extensions. Passing a tree
directly returns a tree-only partial bundle with `log === nothing`.
"""
function with_alignment(rng::AbstractRNG, outbreak::Outbreak, site_model)
    outbreak.tree === nothing &&
        throw(ArgumentError("Cannot simulate an alignment before the tree stage is populated. Use with_tree first."))
    aln = SeqSim.simulate_alignment(rng, outbreak.tree, site_model)
    return Outbreak(outbreak.log, outbreak.tree, aln)
end

with_alignment(outbreak::Outbreak, site_model) = with_alignment(Random.default_rng(), outbreak, site_model)
with_alignment(rng::AbstractRNG, tree, site_model) = with_alignment(rng, Outbreak(nothing, tree, nothing), site_model)
with_alignment(tree, site_model) = with_alignment(Random.default_rng(), tree, site_model)

"""
    outbreak_alignment(rng, tree, site_model)
    outbreak_alignment(tree, site_model)

Return only the tip alignment simulated by SeqSim.
"""
outbreak_alignment(rng::AbstractRNG, tree, site_model) = SeqSim.simulate_alignment(rng, tree, site_model)
outbreak_alignment(tree, site_model) = SeqSim.simulate_alignment(Random.default_rng(), tree, site_model)

"""
    simulate_outbreak_alignment(rng, log, site_model; validate=true)
    simulate_outbreak_alignment(log, site_model; validate=true)

Run the first composed workflow: event log to sampled-ancestry tree to tip
alignment. Returns an `Outbreak` containing the log, tree, and alignment.
"""
function simulate_outbreak_alignment(
    rng::AbstractRNG,
    log,
    site_model;
    validate::Bool=true,
)
    return with_alignment(rng, with_tree(log; validate), site_model)
end

function simulate_outbreak_alignment(log, site_model; validate::Bool=true)
    return simulate_outbreak_alignment(Random.default_rng(), log, site_model; validate)
end

end # module Outbreak
