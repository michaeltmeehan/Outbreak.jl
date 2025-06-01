module Outbreak

using Lazy
using StatsBase

using EpiSim


export AbstractModel, Model
export MTBDModel, BDModel, SIRModel, SEIRModel, SuperSpreaderModel
export MTBDParameters, BDParameters, SIRParameters, SEIRParameters, SuperSpreaderParameters
export AgenticMTBDState, AgenticBDState, AgenticSIRState, AgenticSEIRState, AgenticSuperSpreaderState
export AggregateMTBDState, AggregateBDState, AggregateSIRState, AggregateSEIRState, AggregateSuperSpreaderState

export Seed, Transmission, Recovery, Sampling, Activation
export n_sampled, n_recovered, n_transmissions, n_seeds, n_activations, n_events

export simulate, event_counts

# export BirthDeathModel, MultiTypeBirthDeathModel, SIRModel, SEIRModel, SuperSpreaderModel
# export simulate_events

using SeqSim

export SiteModel, JC, F81, K2P, HKY, GTR, Sequence, Alignment, SequencePropagator
export rand_seq, rand_seq_int

using Random

include("Node.jl")

export AbstractNode, Node, Tree
export isleaf, isbinary, isroot

include("simulate.jl")

export simulate_alignment

include("processing.jl")

# export filter_event_log
export get_sampled_tree, get_tree_stats, get_ltt

include("export.jl")

export build_newick, write_nexus

using EzXML

include("beast/BEASTPlugIns.jl")
using .BEASTPlugIns

export build_elementnode, to_elementnode, make_datetrait_element, make_tree_element, make_beast_xml
export get_tree_stats

end # module Outbreak
