using Outbreak
using Random
using Plots
using StatsPlots
using StatsBase
using DifferentialEquations

λ = 3.0
μ = 0.5
ψ = 0.5
Smax = 150

model = BDModel(; birth_rate=λ, death_rate=μ, sampling_rate=ψ)
rng = Random.MersenneTwister(1234)
ens = simulate(rng, model, 100, stop_condition=s -> s.n_sampled >= Smax)
T = [sim.state_log.t[end] for sim in ens if n_sampled(sim) >= Smax]
Tbar = mean(T)
println("Mean time to reach $Smax samples: $Tbar")
Tpred = 1. / (λ - μ - ψ) * log(1. + (λ - μ - ψ) * Smax / ψ) # Prediction not counting bias due to non-extinction / reaching Smax
println("Predicted time to reach $Smax samples: $Tpred")


function p0(t, λ, μ, ψ, T)
    r = λ - μ - ψ
    c₁ = sqrt(r^2 + 4 * λ * ψ)
    c₂ = -r / c₁
    τ = T - t
    return ((λ + μ + ψ) + c₁ * (exp(-c₁ * τ) * (1. - c₂) - (1. + c₂)) / (exp(-c₁ * τ) * (1 - c₂) + (1 + c₂))) / (2 * λ)
end

# At t = T, p0(T, λ, μ, ψ, T) = 1
# At t → -∞, p0(t, λ, μ, ψ, T) → (λ + μ + ψ - c₁) / (2 * λ)

function dLdt!(dL, L, par, t)
    λ, μ, ψ, T = par
    dL[1] = λ * L[1] * (1. - p0(t, λ, μ, ψ, T)) - ψ * L[1]
end

params = (λ, μ, ψ, Tbar)
prob = ODEProblem(dLdt!, L0, tspan, params)
sol = solve(prob, Tsit5(), saveat=1e-3)


function coupled!(du, u, par, t)
    L, q = u
    λ, μ, ψ, T = par
    du[1] = λ * L * (1. - q) - ψ * L
    du[2] = μ + λ * q^2 + μ - (λ + μ + ψ) * q
end

trees = [get_sampled_tree(sim) for sim in ens if n_sampled(sim) >= 150]

t_max = maximum([tree[1].time for tree in trees])



### --- Predicted lineage-through-time plots --- ###
function ρ(tᵢ, tⱼ, λ, μ, ψ)
    return (μ + ψ - λ) * (tⱼ - tᵢ)
end


function ξ(tᵢ, tⱼ, λ, μ, ψ)
    return (μ + ψ) * (exp(ρ(tᵢ, tⱼ, λ, μ, ψ)) - 1) / ((μ + ψ) * exp(ρ(tᵢ, tⱼ, λ, μ, ψ)) - λ)
end


function η(tᵢ, tⱼ, λ, μ, ψ)
    return λ / (μ + ψ) * ξ(tᵢ, tⱼ, λ, μ, ψ)
end


function p0(tⱼ, tₖ, λ, μ, ψ)
    r = λ - μ - ψ
    c₁ = sqrt(r^2 + 4 * λ * ψ)
    c₂ = -r / c₁
    τ = tₖ - tⱼ
    return ((λ + μ + ψ) + c₁ * (exp(-c₁ * τ) * (1. - c₂) - (1. + c₂)) / (exp(-c₁ * τ) * (1 - c₂) + (1 + c₂))) / (2 * λ)
end

function G0(tⱼ, tₖ, λ, μ, ψ)
    z₁ = (λ + μ + ψ + sqrt((λ + μ + ψ)^2 - 4 * μ * λ)) / (2 * λ)
    z₂ = (λ + μ + ψ - sqrt((λ + μ + ψ)^2 - 4 * μ * λ)) / (2 * λ)
    τ = tₖ - tⱼ
    return (z₁ * (1. - z₂) - z₂ * (1. - z₁) * exp(λ * (z₁ - z₂) * τ)) / ((1. - z₂) - (1. - z₁) * exp(λ * (z₁ - z₂) * τ))
end


function χ(tᵢ, tⱼ, tₖ, λ, μ, ψ)
    ξᵢⱼ = ξ(tᵢ, tⱼ, λ, μ, ψ)
    ηᵢⱼ = η(tᵢ, tⱼ, λ, μ, ψ)
    p₀ⱼₖ = p0(tⱼ, tₖ, λ, μ, ψ)
    return 1. - (1. - p₀ⱼₖ) * (1. - ξᵢⱼ) / (1. - p₀ⱼₖ * ηᵢⱼ)
end


function ϕ(tᵢ, tⱼ, tₖ, λ, μ, ψ)
    ηᵢⱼ = η(tᵢ, tⱼ, λ, μ, ψ)
    p₀ⱼₖ = p0(tⱼ, tₖ, λ, μ, ψ)
    return (1. - p₀ⱼₖ) * ηᵢⱼ / (1 - p₀ⱼₖ * ηᵢⱼ)
end



function pn(n, tᵢ, tⱼ, λ, μ, ψ)
    ξᵢⱼ = ξ(tᵢ, tⱼ, λ, μ, ψ)
    n == 0 && return ξᵢⱼ
    ηᵢⱼ = η(tᵢ, tⱼ, λ, μ, ψ)
    return (1. - ξᵢⱼ) * (1. - ηᵢⱼ) * ηᵢⱼ^(n-1)
end


function ps(s, tᵢ, tⱼ, tₖ, λ, μ, ψ)
    χᵢⱼₖ = χ(tᵢ, tⱼ, tₖ, λ, μ, ψ)
    s == 0 && return χᵢⱼₖ
    ϕᵢⱼₖ = ϕ(tᵢ, tⱼ, tₖ, λ, μ, ψ)
    return (1. - χᵢⱼₖ) * (1. - ϕᵢⱼₖ) * ϕᵢⱼₖ^(s-1)
end


function get_distribution(a::Matrix{Int}, max::Int)
    n = size(a, 1)
    m = size(a, 2)
    dist = zeros(n, max+1)
    for i in 1:n
        freqs = countmap(a[i, :])
        for (k, v) in freqs
            dist[i, k + 1] = v / m
        end
    end
    return dist
end


λ = 3.0
μ = 0.5
ψ = 0.5
t_max = 1.0
t_span = collect(0.:0.25:t_max)
n_sim = 50_000

rng = Random.MersenneTwister(1234)
model = BDModel(; birth_rate=λ, death_rate=μ, sampling_rate=ψ)
ens = simulate(rng, model, n_sim, stop_condition=s -> s.t >= t_max)
trees = get_sampled_tree.(ens)
ltt = reduce(hcat, get_ltt.(trees, Ref(t_span)))
# mean_ltt = mean(ltt, dims=2)[:]
# plot(t_span, mean_ltt)
# predicted_ltt = (1. .- p0.(t_span, λ, μ, ψ, t_max)) .* exp.((λ - μ - ψ) .* t_span)
# [predicted_ltt, mean_ltt]


N = reduce(hcat, [get_state(sim, t_span).I for sim in ens])
max_N = maximum(N)
N_dist = get_distribution(N, max_N)
N_theor = reduce(vcat, [pn.(0:max_N, 0., t, λ, μ, ψ)' for t in t_span])
# N_theor .= round.(N_theor, digits=Int(log10(n_sim)))


max_L = maximum(ltt)
L_dist = get_distribution(ltt, max_L)
L_theor = reduce(vcat, [ps.(0:max_L, 0., t, t_max, λ, μ, ψ)' for t in t_span])
L_theor .= round.(L_theor, digits=Int(log10(n_sim)))



#----- Joint distribution with S -----#
function η_0(tᵢ, tⱼ, λ, μ, ψ)
    b = λ + μ + ψ
    Δ = sqrt(b^2 - 4. * λ * μ)
    r₁ = (b + Δ) / (2. * λ)
    r₂ = (b - Δ) / (2. * λ)
    τ = tⱼ - tᵢ
    return 2* λ * (1. - exp(-Δ * τ)) / ((b + Δ) - (b - Δ) * exp(-Δ * τ))
end

η_0ij = (t) -> η_0(tᵢ, t, λ, μ, ψ)
ts = range(0.0, stop=tₗ, length=100)
dηdt = ForwardDiff.derivative.(η_0ij, ts)
[λ .- (λ + μ + ψ) * η_0ij.(ts) .+ μ * η_0ij.(ts).^2 dηdt]

function ζ_0(tᵢ, tⱼ, λ, μ, ψ)
    b = λ + μ + ψ
    Δ = sqrt(b^2 - 4. * λ * μ)
    r₁ = (b + Δ) / (2. * λ)
    r₂ = (b - Δ) / (2. * λ)
    τ = tⱼ - tᵢ
    return 4 * Δ^2 * exp(-Δ * τ) / ((b + Δ) - (b - Δ) * exp(-Δ * τ))^2
end

ζ_0ij = (t) -> ζ_0(tᵢ, t, λ, μ, ψ)
ts = range(0.0, stop=tₗ, length=100)
dζdt = ForwardDiff.derivative.(ζ_0ij, ts)
[(-(λ + μ + ψ) .+ 2. * μ .* η_0ij.(ts)) .* ζ_0ij.(ts) dζdt]


function π_0(tᵢ, tⱼ, λ, μ, ψ)
    return μ / λ * η_0(tᵢ, tⱼ, λ, μ, ψ)
    # b = λ + μ + ψ
    # Δ = sqrt(b^2 - 4. * λ * μ)
    # r₁ = (b + Δ) / (2. * λ)
    # r₂ = (b - Δ) / (2. * λ)
    # τ = tⱼ - tᵢ
    # return 2 * λ / Δ * log(((b + Δ) - (b - Δ) * exp(-Δ * τ)) / (2 * Δ))
end


π_0ij = (t) -> π_0(tᵢ, t, λ, μ, ψ)
ts = range(0.0, stop=tₗ, length=100)
dπdt = ForwardDiff.derivative.(π_0ij, ts)
[μ .*ζ_0ij.(ts) dπdt]


[π_0ij.(ts) .+ ζ_0ij.(ts) ./ (1. .- η_0ij.(ts)) p0.(tᵢ, ts, λ, μ, ψ)]

function ξ_0(tᵢ, tⱼ, λ, μ, ψ)
    return π_0(tᵢ, tⱼ, λ, μ, ψ) / p0(tᵢ, tⱼ, λ, μ, ψ)
end


function pns0(n, tᵢ, tⱼ, λ, μ, ψ)
    n == 0 && return π_0(tᵢ, tⱼ, λ, μ, ψ)
    return ζ_0(tᵢ, tⱼ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ)^(n-1)
end


function pngs0(n, tᵢ, tⱼ, λ, μ, ψ)
    n == 0 && return ξ_0(tᵢ, tⱼ, λ, μ, ψ)
    return (1. - ξ_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. - η_0(tᵢ, tⱼ, λ, μ, ψ)) * η_0(tᵢ, tⱼ, λ, μ, ψ)^(n-1)
end


S = reduce(hcat, [Outbreak.get_S(sim, t_span) for sim in ens])


P = fill(0.,length(t_span), max_N)
n_0 = fill(0, length(t_span))
for i in eachindex(t_span)
    for j in 1:n_sim
        if S[i, j] == 0
            n_0[i] += 1
            P[i, N[i,j]+1] += 1.
        end
    end
end
P ./ n_sim

P_theor = reduce(vcat, [pns0.(0:max_N, 0., t, λ, μ, ψ)' for t in t_span])

[sum(P ./ n_sim, dims=2) p0.(0., t_span, λ, μ, ψ)]

P ./ n_0

P_theor = reduce(vcat, [pngs0.(0:max_N, 0., t, λ, μ, ψ)' for t in t_span])



# Extract trees with A_i^(jl) = 1
# That is, trees that start with a single individual that goes unsampled until time tⱼ, but does get sampled before tₗ
sims_that_match = fill(0, n_sim)
tj_index = 3
tl_index = 5
for i in 1:n_sim
    if S[tj_index, i] == 0 && S[tl_index, i] > 0
        sims_that_match[i] = 1
    end
end


N_under_condition = Vector{Int}()
ltt_under_condition = Vector{Int}()
for i in 1:n_sim
    if sims_that_match[i] == 1
        push!(N_under_condition, N[tj_index, i])
        push!(ltt_under_condition, ltt[tj_index, i])
    end
end

N_under_condition_dist = fill(0., maximum(N_under_condition) + 1)
ltt_under_condition_dist = fill(0., maximum(ltt_under_condition))
for n in N_under_condition
    N_under_condition_dist[n + 1] += 1.
end
for l in ltt_under_condition
    ltt_under_condition_dist[l] += 1.
end
N_under_condition_dist ./= sum(sims_that_match)
ltt_under_condition_dist ./= sum(sims_that_match)


nx = collect(0:maximum(N_under_condition))
px = p0(t_span[tj_index], t_span[tl_index], λ, μ, ψ)
ηx = η_0(0., t_span[tj_index], λ, μ, ψ)

[(1. .- px.^nx) .* (1. .- ηx) .* ηx.^(nx .- 1) .* (1. .- px .* ηx) ./ (1. .- px) N_under_condition_dist]

lx = collect(1:maximum(ltt_under_condition))
ηy = (1. .- px) .* ηx ./ (1. .- px .* ηx)

[(1. .- ηy) .* ηy.^(lx .- 1) ltt_under_condition_dist]


1. - ξ(tᵢ, tₖ, λ, μ, ψ) ≈ (1. - ξ(tᵢ, tⱼ, λ, μ, ψ)) * (1. - ξ(tⱼ, tₖ, λ, μ, ψ)) / (1. - ξ(tⱼ, tₖ, λ, μ, ψ) * η(tᵢ, tⱼ, λ, μ, ψ))


# (1. - ξ_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. - ξ_0(tⱼ, tₖ, λ, μ, ψ)) / (1. - ξ_0(tᵢ, tₖ, λ, μ, ψ)) ≈ (1. - η_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. - η_0(tⱼ, tₖ, λ, μ, ψ)) / (1. - η_0(tᵢ, tₖ, λ, μ, ψ))

# 1. - ξ_0(tᵢ, tₖ, λ, μ, ψ) ≈ (1. - ξ_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. - ξ_0(tⱼ, tₖ, λ, μ, ψ)) / (1. - ξ_0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ))

# 1. - ξ_0(tᵢ, tₖ, λ, μ, ψ) ≈ (1. - ξ_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. - π_0(tⱼ, tₖ, λ, μ, ψ)) / (1. - π_0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ))

# ξ_0(tᵢ, tₖ, λ, μ, ψ) ≈ ξ_0(tⱼ, tₖ, λ, μ, ψ) * π_0(tᵢ, tⱼ, λ, μ, ψ) +  (ζ_0(tᵢ, tⱼ, λ, μ, ψ) * π_0(tⱼ, tₖ, λ, μ, ψ)) / (1. - π_0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ))

# (1. - ξ_0(tᵢ, tₖ, λ, μ, ψ)) * (1. - η_0(tᵢ, tₖ, λ, μ, ψ)) ≈ (1. - ξ_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. - η_0(tᵢ, tⱼ, λ, μ, ψ)) * p0(tⱼ, tₖ, λ, μ, ψ) / (p0(tᵢ, tₖ, λ, μ, ψ) * (1. - π_0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ))^2)

# (1. - ξ_0(tᵢ, tₖ, λ, μ, ψ)) * (1. - η_0(tᵢ, tₖ, λ, μ, ψ)) ≈ ζ_0(tⱼ, tₖ, λ, μ, ψ) / (p0(tᵢ, tₖ, λ, μ, ψ) * (1. - ξ_0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ))^2)


ζ_0(tᵢ, tⱼ, λ, μ, ψ) ≈ (1. - ξ_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. - η_0(tᵢ, tⱼ, λ, μ, ψ)) * p0(tᵢ, tⱼ, λ, μ, ψ)


ζ_0(tᵢ, tₖ, λ, μ, ψ) ≈ ζ_0(tᵢ, tⱼ, λ, μ, ψ) * ζ_0(tⱼ, tₖ, λ, μ, ψ) / (1. - π_0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ))^2

1. - η_0(tᵢ, tₖ, λ, μ, ψ) ≈ (1. - p0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. - η_0(tⱼ, tₖ, λ, μ, ψ)) / ((1. - π_0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ)))

function α_0(tᵢ, tⱼ, tₖ, λ, μ, ψ)
    p0_jk = p0(tⱼ, tₖ, λ, μ, ψ)
    ξ0_jk = ξ_0(tⱼ, tₖ, λ, μ, ψ)
    η0_ij = η_0(tᵢ, tⱼ, λ, μ, ψ)
    return 1. - p0_jk * (1. - ξ0_jk) * (1. - η0_ij) / (1. - p0_jk * η0_ij)
end

1. - η_0(tᵢ, tₖ, λ, μ, ψ) ≈ (1. - η_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. - η_0(tⱼ, tₖ, λ, μ, ψ)) / (1. - α_0(tᵢ, tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ))

π_0(tᵢ, tₖ, λ, μ, ψ) ≈ π_0(tᵢ, tⱼ, λ, μ, ψ) + ζ_0(tᵢ, tⱼ, λ, μ, ψ) * π_0(tⱼ, tₖ, λ, μ, ψ) / (1. - π_0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ))

p0(tᵢ, tₖ, λ, μ, ψ) ≈ π_0(tᵢ, tⱼ, λ, μ, ψ) + ζ_0(tᵢ, tⱼ, λ, μ, ψ) * p0(tⱼ, tₖ, λ, μ, ψ) / (1. - p0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ))


p0(tᵢ, tⱼ, λ, μ, ψ) - p0(tᵢ, tₖ, λ, μ, ψ) ≈ (p0(tᵢ, tⱼ, λ, μ, ψ) - π_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. - p0(tⱼ, tₖ, λ, μ, ψ)) / (1. - p0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ))


1 - π_0(tᵢ, tₖ, λ, μ, ψ) / p0(tᵢ, tⱼ, λ, μ, ψ) ≈ (1. - ξ_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. - π_0(tⱼ, tₖ, λ, μ, ψ)) / (1. - π_0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ))


function η_0(tᵢ, tⱼ, tₖ, λ, μ, ψ)
    p0_jk = p0(tⱼ, tₖ, λ, μ, ψ)
    η0_ij = η_0(tᵢ, tⱼ, λ, μ, ψ)
    return (1. - p0_jk) * η0_ij / (1. - p0_jk * η0_ij)
end

# 1. - η_0(tᵢ, tⱼ, tₗ, λ, μ, ψ) ≈ (1. - η_0(tᵢ, tₗ, λ, μ, ψ)) / (1. - η_0(tⱼ, tₗ, λ, μ, ψ))


ψ * ζ_0(tᵢ, tⱼ, λ, μ, ψ) / (1. - η_0(tᵢ, tⱼ, λ, μ, ψ))^2 ≈ -μ + (λ + μ + ψ) * p0(tᵢ, tⱼ, λ, μ, ψ) - λ * p0(tᵢ, tⱼ, λ, μ, ψ)^2



# C = μ + ψ * π_0(tᵢ, tⱼ, λ, μ, ψ) / (λ * (1. - η_0(tᵢ, tⱼ, λ, μ, ψ)))
# B = (λ + μ + ψ) + ψ / (λ * (1. - η_0(tᵢ, tⱼ, λ, μ, ψ)))

# p₊ = (B + sqrt(B^2 - 4. * C)) / 2.
# p₋ = (B - sqrt(B^2 - 4. * C)) / 2.

# p0(tᵢ, tⱼ, λ, μ, ψ)

#-------- Conditioning on at least one sampling --------#
numerator = 0
denominator = 0
for i in 1:n_sim
    if S[tj_index, i] > 0
        denominator += 1
        if N[tj_index, i] == 0
            numerator += 1
        end
    end
end
numerator / denominator
(ξ(tᵢ, tⱼ, λ, μ, ψ) - π_0(tᵢ, tⱼ, λ, μ, ψ)) / (1. - p0(tᵢ, tⱼ, λ, μ, ψ))

numerator = 0
denominator = 0
for i in 1:n_sim
    if S[tj_index, i] > 0
        denominator += 1
        if N[tj_index, i] == 2
            numerator += 1
        end
    end
end
numerator / denominator
((1. - ξ(tᵢ, tⱼ, λ, μ, ψ)) * (1. - η(tᵢ, tⱼ, λ, μ, ψ)) * η(tᵢ, tⱼ, λ, μ, ψ) - ζ_0(tᵢ, tⱼ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ)) / (1. - p0(tᵢ, tⱼ, λ, μ, ψ))



p0(tᵢ, tₖ, λ, μ, ψ) - π_0(tᵢ, tₖ, λ, μ, ψ) ≈ (p0(tᵢ, tⱼ, λ, μ, ψ) - π_0(tᵢ, tⱼ, λ, μ, ψ)) * (p0(tⱼ, tₖ, λ, μ, ψ) - π_0(tⱼ, tₖ, λ, μ, ψ)) * (1. - η_0(tᵢ, tⱼ, λ, μ, ψ)) / ((1. - π_0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ)) * (1. -  p0(tⱼ, tₖ, λ, μ, ψ) * η_0(tᵢ, tⱼ, λ, μ, ψ)))


c = π_0(tⱼ, tₖ, λ, μ, ψ) * p0(tⱼ, tₖ, λ, μ, ψ) - (1. - π_0(tⱼ, tₖ, λ, μ, ψ)) * (1. - p0(tⱼ, tₖ, λ, μ, ψ)) / (1. - η_0(tᵢ, tⱼ, λ, μ, ψ))


#----- Conditioning on at least one sampling -----#
function θ(tᵢ, tⱼ, tₖ, λ, μ, ψ)
    ξᵢⱼ = ξ(tᵢ, tⱼ, λ, μ, ψ)
    ηᵢⱼ = η(tᵢ, tⱼ, λ, μ, ψ)
    p₀ᵢₖ = p0(tᵢ, tₖ, λ, μ, ψ)
    p₀ⱼₖ = p0(tⱼ, tₖ, λ, μ, ψ)
    return 1. - (1. - ξᵢⱼ) * (1. -p₀ⱼₖ) / ((1. - p₀ⱼₖ * ηᵢⱼ) * (1. - p₀ᵢₖ))
end


function ps_cond(s, tᵢ, tⱼ, tₖ, λ, μ, ψ)
    θᵢⱼₖ = θ(tᵢ, tⱼ, tₖ, λ, μ, ψ)
    s == 0 && return θᵢⱼₖ
    ϕᵢⱼₖ = ϕ(tᵢ, tⱼ, tₖ, λ, μ, ψ)
    return (1. - θᵢⱼₖ) * (1. - ϕᵢⱼₖ) * ϕᵢⱼₖ^(s-1)
end


conditioned_trees = [get_sampled_tree(sim) for sim in ens if n_sampled(sim) > 0]
conditioned_ltt = reduce(hcat, get_ltt.(conditioned_trees, Ref(t_span)))
max_L_conditioned = maximum(conditioned_ltt)
L_cond_dist = get_distribution(conditioned_ltt, max_L_conditioned)

L_cond_theor = reduce(vcat, [ps_cond.(0:max_L_conditioned, 0., t, t_max, λ, μ, ψ)' for t in t_span])


#----- Conditioning on 1 lineage at time t_j -----#
tᵢ = 0.0
tⱼ = 0.5
tₖ = 0.75
tₗ = t_max
jtrees = 0
ktrees = 0
for i in 1:n_sim
    if ltt[3, i] == 1
        jtrees += 1
        if ltt[4, i] == 0
        ktrees += 1
        end
    end
end
ktrees / jtrees
θ(tⱼ, tₖ, tₗ, λ, μ, ψ)

## Identities ##
1. - χ(tᵢ, tⱼ, tₗ, λ, μ, ψ) == (1. - p0(tⱼ, tₗ, λ, μ, ψ)) * (1. - ξ(tᵢ, tⱼ, λ, μ, ψ)) / (1. - p0(tⱼ, tₗ, λ, μ, ψ) * η(tᵢ, tⱼ, λ, μ, ψ))
1. - ϕ(tᵢ, tⱼ, tₗ, λ, μ, ψ) ≈ (1. - η(tᵢ, tⱼ, λ, μ, ψ)) / (1. - p0(tⱼ, tₗ, λ, μ, ψ) * η(tᵢ, tⱼ, λ, μ, ψ))

ϕ(tᵢ, tⱼ, tₗ, λ, μ, ψ) ≈ (1. - p0(tⱼ, tₗ, λ, μ, ψ)) * η(tᵢ, tⱼ, λ, μ, ψ) / (1. - p0(tⱼ, tₗ, λ, μ, ψ) * η(tᵢ, tⱼ, λ, μ, ψ))

1. - χ(tᵢ, tⱼ, tₗ, λ, μ, ψ) ≈ ϕ(tᵢ, tⱼ, tₗ, λ, μ, ψ) / η(tᵢ, tⱼ, λ, μ, ψ) * (1. - ξ(tᵢ, tⱼ, λ, μ, ψ))

1. - χ(tᵢ, tₖ, tₗ, λ, μ, ψ) == (1. - ξ(tᵢ, tⱼ, λ, μ, ψ)) * (1. - χ(tⱼ, tₖ, tₗ, λ, μ, ψ)) / (1. - χ(tⱼ, tₖ, tₗ, λ, μ, ψ) * η(tᵢ, tⱼ, λ, μ, ψ))
1. - ϕ(tᵢ, tₖ, tₗ, λ, μ, ψ) == (1. - η(tᵢ, tⱼ, λ, μ, ψ)) * (1. - ϕ(tⱼ, tₖ, tₗ, λ, μ, ψ)) / (1. - χ(tⱼ, tₖ, tₗ, λ, μ, ψ) * η(tᵢ, tⱼ, λ, μ, ψ))


1. - χ(tᵢ, tₖ, tₗ, λ, μ, ψ) ≈ (1. - χ(tⱼ, tₖ, tₗ, λ, μ, ψ)) * (1. - χ(tᵢ, tⱼ, tₗ, λ, μ, ψ)) / (1. - χ(tⱼ, tₖ, tₗ, λ, μ, ψ) * η(tᵢ, tⱼ, λ, μ, ψ)) * η(tᵢ, tⱼ, λ, μ, ψ) / ϕ(tᵢ, tⱼ, tₗ, λ, μ, ψ)

(1. - χ(tᵢ, tⱼ, tₗ, λ, μ, ψ)) * (1. - χ(tⱼ, tₖ, tₗ, λ, μ, ψ)) / (1. - χ(tᵢ, tₖ, tₗ, λ, μ, ψ)) ≈ (1. - p0(tⱼ, tₗ, λ, μ, ψ)) * (1. - ϕ(tⱼ, tₖ, tₗ, λ, μ, ψ)) * (1. - ϕ(tᵢ, tⱼ, tₗ, λ, μ, ψ)) / (1. - ϕ(tᵢ, tₖ, tₗ, λ, μ, ψ))


(1. - χ(tᵢ, tⱼ, tₗ, λ, μ, ψ)) / (1. - ϕ(tᵢ, tⱼ, tₗ, λ, μ, ψ)) ≈ (1. - p0(tⱼ, tₗ, λ, μ, ψ)) * (1. - ξ(tᵢ, tⱼ, λ, μ, ψ)) / (1. - η(tᵢ, tⱼ, λ, μ, ψ))

1. - θ(tᵢ, tⱼ, tₗ, λ, μ, ψ) == (1. - χ(tᵢ, tⱼ, tₗ, λ, μ, ψ)) / (1. - p0(tᵢ, tₗ, λ, μ, ψ))

1. - χ(tᵢ, tⱼ, tₗ, λ, μ, ψ) - ϕ(tᵢ, tⱼ, tₗ, λ, μ, ψ) ≈ (1. - p0(tⱼ, tₗ, λ, μ, ψ)) * (1. - ξ(tᵢ, tⱼ, λ, μ, ψ) - η(tᵢ, tⱼ, λ, μ, ψ)) / (1. - p0(tⱼ, tₗ, λ, μ, ψ) * η(tᵢ, tⱼ, λ, μ, ψ))


p0(tᵢ, tₖ, λ, μ, ψ) ≈ p0(tᵢ, tⱼ, λ, μ, ψ) * p0(tⱼ, tₖ, λ, μ, ψ) * (1. - ξ(tᵢ, tⱼ, λ, μ, ψ)) * (1 - η(tᵢ, tⱼ, λ, μ, ψ)) / (1. - p0(tⱼ, tₖ, λ, μ, ψ) * η(tᵢ, tⱼ, λ, μ, ψ))

# Derivative Identities
using ForwardDiff
χᵢⱼₖ = (t) -> χ(tᵢ, t, tₗ, λ, μ, ψ)
ts = range(0.0, stop=tₗ, length=100)
plot(ts, χᵢⱼₖ.(ts), label="χ(tᵢ, t, tₗ)")

dχdt = ForwardDiff.derivative.(χᵢⱼₖ, ts)
plot(ts, dχdt, label="dχ/dt")
ψ ./ (1. .- p0.(ts, tₗ, λ, μ, ψ)) .* (1. .- χ.(tᵢ, ts, tₗ, λ, μ, ψ)) .* (1. .- ϕ.(tᵢ, ts, tₗ, λ, μ, ψ))

ψ ./ (1. .- p0.(ts, tₗ, λ, μ, ψ) .* η.(tᵢ, ts, λ, μ, ψ)) .* (1. .- ξ.(tᵢ, ts, λ, μ, ψ)) .* (1. .- ϕ.(tᵢ, ts, tₗ, λ, μ, ψ))


ϕᵢⱼₖ = (t) -> ϕ(tᵢ, t, tₗ, λ, μ, ψ)
dϕdt = ForwardDiff.derivative.(ϕᵢⱼₖ, ts)
(λ * (1. .- p0.(ts, tₗ, λ, μ, ψ)) .- ψ * η.(tᵢ, ts, λ, μ, ψ) ./ (1. .- p0.(ts, tₗ, λ, μ, ψ) .* η.(tᵢ, ts, λ, μ, ψ))) .* (1. .- ϕ.(tᵢ, ts, tₗ, λ, μ, ψ))
(λ * (1. .- p0.(ts, tₗ, λ, μ, ψ)) .- ψ * ϕ.(tᵢ, ts, tₗ, λ, μ, ψ) ./ (1. .- p0.(ts, tₗ, λ, μ, ψ))) .* (1. .- ϕ.(tᵢ, ts, tₗ, λ, μ, ψ))

θᵢⱼₖ = (t) -> θ(tᵢ, t, tₗ, λ, μ, ψ)
dθdt = ForwardDiff.derivative.(θᵢⱼₖ, ts)
ψ ./ (1. .- p0.(ts, tₗ, λ, μ, ψ)) .* (1. .- θ.(tᵢ, ts, tₗ, λ, μ, ψ)) .* (1. .- ϕ.(tᵢ, ts, tₗ, λ, μ, ψ))


p1 = (t) -> (1. - θ(tᵢ, t, tₗ, λ, μ, ψ)) * (1. - ϕ(tᵢ, t, tₗ, λ, μ, ψ))
dp1dt = [ForwardDiff.derivative(p1, t) for t in ts]
-(λ * (1. .- p0.(ts, tₗ, λ, μ, ψ)) .- ψ ./ (1. .- p0.(ts, tₗ, λ, μ, ψ)) .* (2. .* ϕ.(tᵢ, ts, tₗ, λ, μ, ψ) .- 1.)) .* (1. .- θ.(tᵢ, ts, tₗ, λ, μ, ψ)) .* (1. .- ϕ.(tᵢ, ts, tₗ, λ, μ, ψ))


#----- Conditioning on survival to time T -----#
function p_cond(n, tᵢ, tⱼ, tₖ, λ, μ, ψ)
    ξᵢⱼ = ξ(tᵢ, tⱼ, λ, μ, ψ)
    ξᵢₖ = ξ(tᵢ, tₖ, λ, μ, ψ)
    ξⱼₖ = ξ(tⱼ, tₖ, λ, μ, ψ)
    ηᵢⱼ = η(tᵢ, tⱼ, λ, μ, ψ)
    n == 0 && return 0.
    return (1. - ξᵢⱼ) * (1. - ηᵢⱼ) / (1. - ξᵢₖ) * (1. - ξⱼₖ^n) * ηᵢⱼ^(n-1)
end


function ps_cond(s, tᵢ, tⱼ, tₖ, λ, μ, ψ)
    ξᵢⱼ = ξ(tᵢ, tⱼ, λ, μ, ψ)
    ξᵢₖ = ξ(tᵢ, tₖ, λ, μ, ψ)
    ξⱼₖ = ξ(tⱼ, tₖ, λ, μ, ψ)
    ηᵢⱼ = η(tᵢ, tⱼ, λ, μ, ψ)
    p₀ⱼₖ = p0(tⱼ, λ, μ, ψ, tₖ)
    s == 0 && return (1. - ηᵢⱼ * ξⱼₖ) * (1. - ηᵢⱼ) * p₀ⱼₖ / ((1. - p₀ⱼₖ * ηᵢⱼ) * (1. - p₀ⱼₖ * ξⱼₖ * ηᵢⱼ))
    return (1. - ξᵢⱼ) * (1. - ηᵢⱼ) / (1. - ξᵢₖ) * (1. - p₀ⱼₖ)^s * ηᵢⱼ^(s-1) * (1. / (1. - p₀ⱼₖ * ηᵢⱼ)^(s+1) - ξⱼₖ^s / (1. - p₀ⱼₖ * ξⱼₖ * ηᵢⱼ)^(s+1))
end


N_cond_theor = reduce(vcat, [p_cond.(0:max_N, 0., t, t_max, λ, μ, ψ)' for t in t_span])
N_cond_theor .= round.(N_cond_theor, digits=Int(log10(n_sim)))

N_conditioned = reduce(hcat, [get_state(sim, t_span).I for sim in ens if sim.state_log.t[end] >= t_max])
max_N_conditioned = maximum(N_conditioned)
N_dist_conditioned = get_distribution(N_conditioned, max_N_conditioned)

conditioned_trees = [get_sampled_tree(sim) for sim in ens if sim.state_log.t[end] >= t_max]
conditioned_ltt = reduce(hcat, get_ltt.(conditioned_trees, Ref(collect(t_span))))

max_L_conditioned = maximum(conditioned_ltt)
L_cond_dist = get_distribution(conditioned_ltt, max_L)

L_cond_theor = reduce(vcat, [ps_cond.(0:max_L, 0., t, t_max, λ, μ, ψ)' for t in t_span])