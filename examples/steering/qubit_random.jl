#In this example, we study the visibility for which random states admit a local-hidden-state model
#We've beforehand generated 100 random steerable two-qubit states using random_state from Ket.jl and critical_radius
#So, we first load such states
using Serialization
states = deserialize("random_steerable_states.dat")

#Now, we compute their lhs visibility, first using the method from [arXiv:1808.09349](https://arxiv.org/abs/1808.09349)
using LpJm
polytope162 = deserialize("../../polytopes/dim2/polytope162cov.dat")
shr162 = shrinking_factor(polytope162)
polytope162_outer = bloch_op.(bloch_vec.(polytope162) / shr162)

lower_162 = critical_radius(states[1], polytope162)
upper_162 = critical_radius(states[1], polytope162_outer)

#Then, we use our methods to compute bounds to the same visibility
measurements = deserialize("../../polytopes/dim2/measurements406.dat")

polytope614 = deserialize("../../polytopes/dim2/polytope614.dat")
shr614 = shrinking_factor(polytope614)
polytope614_outer = bloch_op.(bloch_vec.(polytope614) / shr614)

#The following steps are not mandatory, but increase performance (the same thing is done in critical_radius)
using LinearAlgebra
using Ket
F = kron(I(2), sqrt(inv(partial_trace(states[1], 1, [2, 2]))))
new_state = F * states[1] * F
new_state ./= tr(new_state)

lower = lhs_visibility_lower(new_state, measurements, polytope614)
upper = lhs_visibility_upper(new_state, measurements, polytope614_outer)
