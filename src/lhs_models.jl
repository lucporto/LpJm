function lhs_visibility_upper(
    ρ::AbstractMatrix,
    all_Aax::Vector{<:Measurement},
    all_χ::Vector{<:AbstractMatrix},
    dims::Vector{<:Integer} = Ket._equal_sizes(ρ);
    return_probs::Bool = false,
    solver = Mosek.Optimizer
)
    m = length(all_Aax)
    o = length(all_Aax[1])

    n = length(all_χ)

    model = Model(solver)
    set_silent(model)

    p = [@variable(model, [1:o, 1:m], lower_bound = 0.0) for _ ∈ 1:n]
    @variable(model, η)
    @variable(model, π[1:n] .>= 0)

    noisy_ρ = @expression(model, η * ρ + (1 - η) * kron(I(dims[1]) / dims[1], partial_trace(ρ, 1, dims)))
    for x ∈ 1:m
        @constraint(model, [i in 1:n], sum(p[i][a, x] for a ∈ 1:o) == π[i])
        @constraint(
            model,
            [a in 1:o],
            partial_trace(kron(all_Aax[x][a], I(dims[2])) * noisy_ρ, 1, dims) .==
            sum(p[i][a, x] * Matrix(all_χ[i]) for i ∈ 1:n)
        )
    end
    @constraint(model, η ≤ 2)

    @objective(model, Max, η)

    optimize!(model)

    primal_status(model) == FEASIBLE_POINT || throw(error(raw_status(model)))
    is_solved_and_feasible(model) || @warn("Something went wrong, but the answer is feasible: $(raw_status(model))")

    return return_probs ? (objective_value(model), broadcast(x -> value.(x), p)) : objective_value(model)
end
export lhs_visibility_upper

function lhs_visibility_lower(
    ρ::AbstractMatrix,
    all_Aax::Vector{<:Measurement},
    all_ρ::Vector{<:AbstractMatrix},
    dims::Vector{<:Integer} = Ket._equal_sizes(ρ);
    r = shrinking_factor(all_Aax),
    return_probs::Bool = false,
    solver = Mosek.Optimizer
)
    m = length(all_Aax)
    o = length(all_Aax[1])
    depolarized_Aax =
        [[(1 / r) * all_Aax[x][a] + (1 - 1 / r) * tr(all_Aax[x][a]) * I(dims[1]) / dims[1] for a ∈ 1:o] for x ∈ 1:m]

    return lhs_visibility_upper(ρ, depolarized_Aax, all_ρ; return_probs, solver)
end
export lhs_visibility_lower
