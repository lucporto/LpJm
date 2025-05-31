function lp_steering_robustness(
    σ::Matrix{<:AbstractMatrix},
    all_ρ::Vector{<:AbstractMatrix};
    return_probs::Bool = false,
    solver = Mosek.Optimizer
)
    o, m = size(σ)
    d = size(σ[1, 1], 1)

    n = length(all_ρ)

    model = Model(solver)
    set_silent(model)

    p = [@variable(model, [1:o, 1:m], lower_bound = 0.0) for _ ∈ 1:n]
    @variable(model, η)
    @variable(model, π[1:n] .>= 0)

    for x ∈ 1:m
        @constraint(model, [i in 1:n], sum(p[i][a, x] for a ∈ 1:o) == π[i])
        @constraint(
            model,
            [a in 1:o],
            η * σ[a, x] + (1 - η) * tr(σ[a, x]) * I(d) / d .== sum(p[i][a, x] * Matrix(all_ρ[i]) for i ∈ 1:n)
        )
    end
    @constraint(model, η ≤ 2)

    @objective(model, Max, η)

    optimize!(model)

    primal_status(model) == FEASIBLE_POINT || throw(error(raw_status(model)))
    is_solved_and_feasible(model) || @warn("Something went wrong, but the answer is feasible: $(raw_status(model))")

    return return_probs ? (objective_value(model), broadcast(x -> value.(x), p)) : objective_value(model)
end
export lp_steering_robustness

function lp_incompatibility_robustness(
    Aax::Vector{<:Vector{<:AbstractMatrix}},
    all_ρ::Vector{<:AbstractMatrix};
    return_probs::Bool = false,
    solver = Mosek.Optimizer
)
    σ = [transpose(Aax[x][a]) / size(Aax[1][1], 1) for a ∈ eachindex(Aax[1]), x ∈ eachindex(Aax)]
    return lp_steering_robustness(σ, all_ρ; return_probs, solver)
end
export lp_incompatibility_robustness
