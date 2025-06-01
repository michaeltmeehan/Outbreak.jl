## ---- Birth-death tree statistics ---- ##
using DataFrames
using Random
using Base.Threads
using Outbreak
using AlgebraOfGraphics
using CairoMakie

rng = Random.MersenneTwister(1234)


function simulate_tree_stats(rng::AbstractRNG, model::AbstractModel, n::Int, Smax::Int)
    ens = simulate(rng, model, n, stop_condition= s -> s.n_sampled >= Smax)
    out = DataFrame()
    for sim in ens
        if n_sampled(sim) < Smax
            continue
        end
        tree = get_sampled_tree(sim)
        push!(out, get_tree_stats(tree))
    end
    return out
end


R0_grid = [1.5, 2., 2.5, 3., 3.5, 4.]
δ_grid = [0.1, 0.2, 0.3, 0.4, 0.5]
ψ_grid = [0.1, 0.2, 0.3, 0.4, 0.5]


function run_grid(rng, model_type::Symbol, R0s::Vector{Float64}, δs::Vector{Float64}, ψs::Vector{Float64}; n::Int = 1000, Smax::Int = 50)
    out = DataFrame()
    for (R0, δ, ψ) in Iterators.product(R0s, δs, ψs)
        if model_type == :bd
            model = BDModel(; birth_rate=R0 * (δ + ψ), death_rate=δ, sampling_rate=ψ)
        elseif model_type == :ss
            model = SuperSpreaderModel(; R0=R0, recovery_rate=fill(δ, 2), sampling_rate=fill(ψ, 2))
        else
            error("Unknown model type: $model_type")
        end
        df = simulate_tree_stats(rng, model, n, Smax)
        df.R0 .= R0
        df.δ .= δ
        df.ψ .= ψ
        df.n_sampled .= Smax
        df.model .= model_type
        append!(out, df)
    end
    return out
end

# ---- Super-spreader tree statistics ---- #

out_bd = run_grid(rng, :bd, R0_grid, δ_grid, ψ_grid)
out_ss = run_grid(rng, :ss, R0_grid, δ_grid, ψ_grid)

out = vcat(out_bd, out_ss)


function plot_tree_stat(out::DataFrame, stat::Symbol, label=String)
    data(out) *
        mapping(:R0 => nonnumeric,
                stat .=> label,
                color = :model,
                dodge = :model,
                col = :δ => nonnumeric,
                row = :ψ => nonnumeric
                ) *
        visual(BoxPlot; width=1, markersize=8.) |>
        draw(; figure=(size = (1000, 1000),))
end

plot_tree_stat(out, :tree_height, "Tree height (years)")
plot_tree_stat(out, :colless_index, "Colless index")
plot_tree_stat(out, :sackin_index, "Sackin index")
plot_tree_stat(out, :n_cherries, "Number of cherries")
plot_tree_stat(out, :n_ladder_nodes, "Number of ladder nodes")
plot_tree_stat(out, :il_portion, "Internal ladder portion")
plot_tree_stat(out, :max_width_on_depth, "Max width on depth")
plot_tree_stat(out, :max_width_difference, "Max width difference")
plot_tree_stat(out, :ladder_length, "Ladder length")

data(out) *
    mapping(:R0 => nonnumeric,
            (:sackin_index, :tree_height) => (/) => "(Normalized) Sackin index",
            color = :model,
            dodge = :model,
            col = :δ => nonnumeric => "Recovery rate",
            row = :ψ => nonnumeric => "Sampling rate"
            ) *
    visual(BoxPlot; width=1, markersize=8.) |>
    draw(; figure=(size = (1000, 1000),))



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