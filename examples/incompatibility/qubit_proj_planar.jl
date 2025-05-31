#In this example, we test a conjecture by [arXiv:2003.00785](https://arxiv.org/abs/2003.00785)
using LinearAlgebra
using Ket
using LpJm

function conjecture_lower_bound(α::Vector{<:Real})
    N = length(α)
    α = sort(α)
    return 1 / (sum(sin((α[k] - α[k-1]) / 2) for k ∈ 2:N) + cos(α[N] / 2))
end

function planar_povm(α::Real)
    P = Hermitian(cos(α) * pauli(1) + sin(α) * pauli(2))
    return [I(2) / 2 + P / 2, I(2) / 2 - P / 2]
end

function states_regular_polygon(n::Integer)
    return [cleanup!(bloch_op([cos(α), sin(α), 0])) for α ∈ range(0; step = 2 * pi / n, length = n)]
end

#The idea is to sample a random set of qubit planar projective measurements and try to falsify the conjecture
#If our_bound > conj, the conjecture is false
angles = pi .* rand(100)
angles[1] = 0
conj = conjecture_lower_bound(angles)
our_bound = lp_incompatibility_robustness(planar_povm.(angles), states_regular_polygon(400))
