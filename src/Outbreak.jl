module Outbreak

using EpiSim

export BirthDeathModel, MultiTypeBirthDeathModel, SIRModel, SEIRModel, SuperSpreaderModel
export simulate_outbreak

using SeqSim

export SiteModel, JC, F81, K2P, HKY, GTR, Sequence, SequencePropagator
export rand_seq

using Random

include("Node.jl")

export AbstractNode

include("simulate.jl")

export simulate_alignment

include("processing.jl")

export filter_event_log, get_sampled_tree

include("export.jl")

export build_newick

end # module Outbreak
