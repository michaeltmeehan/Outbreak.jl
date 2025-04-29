module Outbreak

using EpiSim

export BirthDeathModel, MultiTypeBirthDeathModel, SIRModel, SEIRModel
export simulate_outbreak

using SeqSim

export SiteModel, JC, F81, K2P, HKY, GTR, Sequence
export rand_seq

using Random

include("simulate.jl")

export simulate_alignment

include("processing.jl")

export filter_event_log, get_sampled_tree

end # module Outbreak
