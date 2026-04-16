# Outbreak.jl

Outbreak.jl is the thin orchestration layer for the outbreak modelling package
ecosystem. It keeps the underlying scientific objects intact and composes the
settled bridge APIs from the lower-level packages:

- `EpiSim.jl` produces epidemic event logs.
- `TreeSim.jl` owns canonical tree representation and event-log extraction.
- `SeqSim.jl` owns sequence and alignment simulation.

The central object is a small workflow bundle:

```julia
Outbreak(log)
Outbreak(log, tree)
Outbreak(log, tree, aln)
```

Missing stages are represented by `nothing`. The fields hold the original
objects directly:

```julia
out.log
out.tree
out.aln
```

The bundle is a partial workflow state, not a claim that every stage was
produced by Outbreak.jl. For example, tree-first workflows may use
`Outbreak(nothing, tree, aln)` when an existing tree is the starting point.
Use `iscomplete(out)` only when a workflow requires log, tree, and alignment to
all be present.

## Orchestrated Workflows

Populate a sampled-ancestry tree from an event log:

```julia
using Outbreak: Outbreak, with_tree, simulate_outbreak_tree

out = Outbreak(log)
out = with_tree(out)

tree_out = simulate_outbreak_tree(log)
```

Simulate a tip alignment from a populated tree:

```julia
using Outbreak: with_alignment

out = with_alignment(rng, out, site_model)
```

Or start from an existing tree:

```julia
tree_out = with_alignment(rng, tree, site_model)
tree_out.log === nothing
```

Run the first composed path in one call:

```julia
using Outbreak: simulate_outbreak_alignment

out = simulate_outbreak_alignment(rng, log, site_model)
```

These functions delegate to `TreeSim.tree_from_eventlog` and
`SeqSim.simulate_alignment`. Outbreak.jl does not redefine event-log-to-tree or
tree-to-alignment semantics. Stage transitions fail early when the required
upstream stage is missing.

In this recovery workspace, Outbreak.jl uses local sibling package sources for
TreeSim and SeqSim, and its tests load EpiSim so TreeSim's EpiSim extension is
available. Those package extensions are the integration contracts; Outbreak.jl
only composes them.

## Public Surface

The first orchestration surface is intentionally small:

- `Outbreak`
- `has_log`, `has_tree`, `has_alignment`, `iscomplete`
- `with_tree`, `outbreak_tree`
- `with_alignment`, `outbreak_alignment`
- `simulate_outbreak_tree`
- `simulate_outbreak_alignment`

BDUtils integration, inference workflows, calibration, and richer result
hierarchies are deferred until the underlying interfaces need orchestration at
this level.
