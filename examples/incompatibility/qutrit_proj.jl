#In this example, we will study sets of measurements generated according to the fib_qutrit function.
#Construct the measurement set
using LpJm
m = 10
Max = fibonacci_measurements_qutrit(m)

#SDP solution (do not use for more m > 12)
using Ket
using MosekTools
incomp_rob = 1 / (1 + incompatibility_robustness(Max; noise = "depolarizing", solver = Mosek.Optimizer))

#For the LP, choose a polytope...
using Serialization
#using DoubleFloats
#polytope273 = deserialize("../../polytopes/dim3/polytope273.dat") 
#shr273 = 0.831720

polytope751 = deserialize("../../polytopes/dim3/polytope751.dat")
shr751 = 0.883329

#polytope2191 = deserialize("../../polytopes/dim3/polytope2191.dat")
#shr2191 = 0.9066623

#...and compute the approximate incompatibility robustness:
lower_bound = lp_incompatibility_robustness(Max, polytope751)
upper_bound = lower_bound / shr751
