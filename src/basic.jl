function bloch_op(v::AbstractVector{<:Real})
    d = isqrt(length(v) + 1)
    gm = gellmann(d)
    return Hermitian(gm[1] / d + sum(v[i] * gm[i+1] / 2 for i ∈ eachindex(v)))
end
export bloch_op

function bloch_vec(A::AbstractMatrix)
    !(A ≈ A') && throw(ArgumentError("A must be hermitian"))
    !(tr(A) ≈ 1) && throw(ArgumentError("A must have unit trace"))
    d = size(A)[1]
    gm = gellmann(d)
    return [Real(tr(gm[i] * A)) for i ∈ 2:d^2]
end
export bloch_vec

function dichotomic_povm(A::AbstractMatrix)
    return [A, I(size(A, 1)) - A]
end
export dichotomic_povm

#Only works for qubits
function shrinking_factor(vertices::Vector{<:Vector{<:Real}})
    return Polyhedra.maximum_radius_with_center(Polyhedra.doubledescription(Polyhedra.vrep(vertices)), [0, 0, 0])
end
shrinking_factor(all_ρ::Vector{<:AbstractMatrix}) = shrinking_factor(bloch_vec.(all_ρ))
shrinking_factor(all_Aax::Vector{<:Measurement}) = shrinking_factor(reduce(vcat, all_Aax))
export shrinking_factor

function fibonacci_measurements_qubit(m::Integer)
    Max = Vector{Measurement{ComplexF64}}(undef, m)
    ϕ = pi * (sqrt(5) - 1)
    for i ∈ 1:m
        z = 1 - (i - 1) / (m - 1)
        radius = sqrt(1 - z^2)

        θ = ϕ * (i - 1)
        x = radius * cos(θ)
        y = radius * sin(θ)

        Max[i] = dichotomic_povm(bloch_op([x, y, z]))
    end
    return Max
end
export fibonacci_measurements_qubit

function fibonacci_measurements_qutrit(m::Integer)
    Max = Vector{Measurement{ComplexF64}}(undef, m)
    for i ∈ 1:m
        z = 1 - (i - 1) / (m - 1)
        θ = pi * (sqrt(5) - 1) * (i - 1)
        ϕ = pi * (i - 1) / (m - 1)

        M1 = ketbra([sqrt(1 - z^2) * cos(θ), sqrt(1 - z^2) * sin(θ), z * exp(ϕ * im)])
        M2 = ketbra([sin(θ), -cos(θ), 0])
        M3 = I(3) - (M1 + M2)

        Max[i] = [M1, M2, M3]
    end
    return Max
end
export fibonacci_measurements_qutrit

function find_extremal_povm(E::Vector{<:AbstractMatrix}; solver = Mosek.Optimizer)
    P = copy(E)
    is_extremal, Q = find_direction(P; solver)
    step = 1
    while is_extremal == false
        u = maximum_weight(P, Q, step)
        P += u * Q
        is_extremal, Q = find_direction(P; solver)
        step += 1
    end
    return P
end
export find_extremal_povm

function find_direction(E::Vector{<:AbstractMatrix{T}}; solver = Mosek.Optimizer, tol = 10^-10) where {T<:Number}
    d = size(E[1], 1)
    null_outcomes = findall(x -> maximum(eigvals(x)) ≤ tol, E)
    non_null_outcomes = findall(x -> !(maximum(eigvals(x)) ≤ tol), E)
    P = E[non_null_outcomes]
    o = length(P)
    V = ketbra_nonzero_eigvecs.(P; tol)
    N = [size(V[i], 1) for i ∈ 1:o]
    Q = Vector{Matrix{Complex{Ket._solver_type(T)}}}(undef, length(E))
    for complex ∈ (true, false)
        for a ∈ 1:o
            for j ∈ 1:N[a], k ∈ 1:N[a]
                x, Q[non_null_outcomes] = border_direction(V, [a, j, k], complex; solver)
                Q[null_outcomes] = [zeros((d, d)) for i ∈ 1:length(null_outcomes)]
                if x > tol
                    Q = [Q[i] + Q[i]' for i ∈ eachindex(Q)]
                    return false, Q
                end
            end
        end
    end
    return true, Q
end

function border_direction(
    V::Vector{<:AbstractMatrix},
    index::Vector{<:Integer},
    complex::Bool;
    solver = Mosek.Optimizer
)
    o = length(V)
    N = [size(V[i], 1) for i ∈ 1:o]

    model = Model(solver)
    set_silent(model)

    D = [[@variable(model, set = ComplexPlane()) for _ ∈ 1:N[a], _ ∈ 1:N[a]] for a ∈ 1:o]
    Q = [@expression(model, sum(D[a][m, n] * V[a][m, n] for m ∈ 1:N[a], n ∈ 1:N[a])) for a ∈ 1:o]

    @constraint(model, sum(Q) .== 0)
    for a ∈ 1:o
        @constraint(model, real.(D[a]) .≤ 1)
        @constraint(model, imag.(D[a]) .≤ 1)
    end

    complex ? @objective(model, Max, imag(D[index[1]][index[2], index[3]])) :
    @objective(model, Max, real(D[index[1]][index[2], index[3]]))

    optimize!(model)
    is_solved_and_feasible(model) || throw(error(raw_status(model)))

    return objective_value(model), broadcast(x -> value.(x), Q)
end

function maximum_weight(E::Vector{<:AbstractMatrix}, Q::Vector{<:AbstractMatrix}, step::Integer; tol = 10^-16)
    f(u) = sum(sort(vcat(eigvals.(E + u * Q)...))[1:step]) + step * tol
    return find_zero(f, 0.0)
end

function nonzero_eigvecs(A::AbstractMatrix; tol = 10^-10)
    vals, vecs = eigen(A)
    nonzero_vals = findall(x -> abs(x) > tol, vals)
    V = Vector{Vector{eltype(vecs)}}(undef, length(nonzero_vals))
    for i ∈ eachindex(nonzero_vals)
        V[i] = vecs[:, nonzero_vals[i]]
    end
    return V
end

function ketbra_nonzero_eigvecs(A::AbstractMatrix; tol = 10^-10)
    V = nonzero_eigvecs(A; tol)
    n = length(V)
    B = Matrix{Matrix{eltype(V[1])}}(undef, n, n)
    for i ∈ 1:n
        for j ∈ 1:n
            B[i, j] = V[i] * V[j]'
        end
    end
    return B
end
