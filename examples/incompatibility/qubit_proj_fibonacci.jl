#In this example, we will study sets of measurements generated according to the fibonacci_measurements_qubit function
#First, we construct the measurement set
using LpJm
m = 10
Max = fibonacci_measurements_qubit(m)

#SDP solution (do not use for more m > 20)
using Ket
using MosekTools
incomp_rob = 1 / (1 + incompatibility_robustness(Max; noise = "depolarizing", solver = Mosek.Optimizer))

#For the LP, choose a polytope...
using Serialization
polytope614 = deserialize("../../polytopes/dim2/polytope614.dat")
shr614 = 0.99599557

#polytope5534 = deserialize("../../polytopes/dim2/polytope5534.dat") 
#shr5534 = 0.9995537

#polytope8194 = deserialize("../../polytopes/dim2/polytope8194.dat") 
#shr8194 = 0.9996987

#...and compute the approximate incompatibility robustness:
lower_bound = lp_incompatibility_robustness(Max, polytope614)
upper_bound = lower_bound / shr614
