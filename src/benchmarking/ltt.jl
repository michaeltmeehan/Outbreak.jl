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
