function critical_radius(ρ::AbstractMatrix, polytope::Vector{<:Vector{<:Real}}; solver = Mosek.Optimizer)
    normals, offsets = all_planes(polytope)

    F = kron(I(2), sqrt(inv(partial_trace(ρ, 1, [2, 2]))))
    ρ = F * ρ * F
    ρ ./= tr(ρ)

    a = Vector{Float64}(undef, 3)
    T = Matrix{Float64}(undef, 3, 3)
    for i ∈ 1:3
        a[i] = real(tr(pauli(i) * partial_trace(ρ, 2, [2, 2])))
        for j ∈ 1:3
            T[i, j] = real(tr(ρ * kron(pauli(i), pauli(j))))
        end
    end

    model = Model(solver)
    set_silent(model)

    r = @variable(model)
    @variable(model, probs[eachindex(polytope)] .>= 0)

    b = [
        @expression(
            model,
            sum(
                probs[j] * abs(-offsets[i] + normals[i]' * polytope[j]) / norm(-offsets[i] * a + T * normals[i]) for
                j ∈ eachindex(polytope)
            )
        ) for i ∈ eachindex(normals)
    ]
    @constraint(model, [i in eachindex(normals)], b[i] ≥ r)
    @constraint(model, sum(probs) == 1)
    @constraint(model, sum(probs[j] * polytope[j] for j ∈ eachindex(polytope)) .== 0)

    @objective(model, Max, r)

    optimize!(model)

    is_solved_and_feasible(model) || throw(error(raw_status(model)))

    return objective_value(model)
end

critical_radius(ρ::AbstractMatrix, polytope::Vector{<:AbstractMatrix}) = critical_radius(ρ, bloch_vec.(polytope))
export critical_radius

#Given three points in 3d, returns the plane that passes through them
function plane(points::Vector{<:Vector{<:Real}})
    a = points[2] - points[1]
    b = points[3] - points[1]
    normal = cross(a, b)
    offset = dot(points[1], normal)
    return normal, offset
end

#Given a finite set of points, returns all triples composed of such points
function all_triples(points::AbstractVector{T}) where {T}
    n = length(points)
    triples = Vector{Vector{T}}(undef, binomial(n, 3))
    count = 1
    for i ∈ 1:(n-2)
        for j ∈ (i+1):(n-1)
            for k ∈ (j+1):n
                triples[count] = [points[i], points[j], points[k]]
                count += 1
            end
        end
    end
    return triples
end

#Given a finite set of points in 3d, all_planes gets all the planes that pass through at least three of those points
function all_planes(points::Vector{Vector{T}}) where {T<:Real}
    triples = all_triples(points)
    n = length(triples)
    normals = Vector{Vector{T}}(undef, n)
    offsets = Vector{T}(undef, n)
    for i ∈ 1:n
        normals[i], offsets[i] = plane(triples[i])
    end
    return normals, offsets
end
