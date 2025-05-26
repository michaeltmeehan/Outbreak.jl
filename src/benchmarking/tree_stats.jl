ss_model = SuperSpreaderModel(; R0=3.)
ss_ens = simulate(ss_model, 1_000, stop_condition = (state) -> state.n_sampled >= 100)
ss_tree_heights = Vector{Float64}(); ss_colless_indices = Vector{Int64}(); ss_sackin_indices = Vector{Float64}();

for sim in ss_ens                                                                                                                                                                          
    events = event_counts(sim)                                                                                                                                                              
    if haskey(events, Sampling) && events[Sampling] == 100                                                                                                                                  
        tree = get_sampled_tree(sim)                                                                                                                                                            
        tree_height, colless_index, sackin_index, _, _, _, _ = get_tree_stats(tree)                                                                                                             
        push!(ss_tree_heights, tree_height)                                                                                                                                                        
        push!(ss_colless_indices, colless_index)                                                                                                                                                   
        push!(ss_sackin_indices, sackin_index)                                                                                                                                                     
    end                                                                                                                                                                                     
end


bd_model = BDModel(; birth_rate=3.0, death_rate=0.9, sampling_rate=0.1)
bd_ens = simulate(bd_model, 1_000, stop_condition = (state) -> state.n_sampled >= 100)
bd_tree_heights = Vector{Float64}(); bd_colless_indices = Vector{Int64}(); bd_sackin_indices = Vector{Float64}();

for sim in bd_ens                                                                                                                                                                          
    events = event_counts(sim)                                                                                                                                                              
    if haskey(events, Sampling) && events[Sampling] == 100                                                                                                                                  
        tree = get_sampled_tree(sim)                                                                                                                                                            
        tree_height, colless_index, sackin_index, _, _, _, _ = get_tree_stats(tree)                                                                                                             
        push!(bd_tree_heights, tree_height)                                                                                                                                                        
        push!(bd_colless_indices, colless_index)                                                                                                                                                   
        push!(bd_sackin_indices, sackin_index)                                                                                                                                                     
    end                                                                                                                                                                                     
end

using StatsPlots
heights = vcat(ss_tree_heights, bd_tree_heights)
colless_indices = vcat(ss_colless_indices, bd_colless_indices)
sackin_indices = vcat(ss_sackin_indices, bd_sackin_indices)
models = vcat(fill("ss", length(ss_tree_heights)), fill("bd", length(bd_tree_heights)))

boxplot(models, heights, legend=true)
boxplot(models, colless_indices, legend=true)
boxplot(models, sackin_indices, legend=true)