#In this example, we study the incompatibility robustness of random extremal qubit POVMs and random qubit projective measurements
#First, we generate the random measurements
using LpJm
using Ket
m = 100
rand_povms = [find_extremal_povm(random_povm(2, 4)) for _ ∈ 1:m]
rand_proj = [find_extremal_povm(random_povm(2, 2)) for _ ∈ 1:m]

#Then, we compute their incompatibility robustness with the LPs
using Serialization
polytope614 = deserialize("../../polytopes/dim2/polytope614.dat")
shr614 = 0.99599557

lower_povm = lp_incompatibility_robustness(rand_povms, polytope614)
upper_povm = lower_povm / shr614

lower_proj = lp_incompatibility_robustness(rand_proj, polytope614)
upper_proj = lower_proj / shr614
