### Code to accompany: *[Measurement Incompatibility and Quantum Steering via Linear Programming](http://arxiv.org/abs/2506.03045)*
#### Lucas Porto, Sébastien Designolle, Sebastian Pokutta and Marco Túlio Quintino

This repository contains an implementation of the methods proposed in the article, which was used to study the examples discussed therein.

The main functions of this repository are:

 - [lp_steering_robustness](https://github.com/lucporto/LpJm/blob/main/src/lp_robustness.jl): Computes bounds on the depolarizing steering robustness of an input assemblage, using a polytope approximation of the set of quantum states that has to given as the second argument of the function.
 - [lp_incompatibility_robustness](https://github.com/lucporto/LpJm/blob/main/src/lp_robustness.jl): Computes bounds on the depolarizing incompatibility robustness of an input measurement set, using a polytope approximation of the set of quantum states that has to given as the second argument of the function.
 - [lhs_visibility_upper](https://github.com/lucporto/LpJm/blob/main/src/lhs_models.jl): Computes an upper bound on the visibility for which a state becomes unsteerable.
 - [lhs_visibility_lower](https://github.com/lucporto/LpJm/blob/main/src/lhs_models.jl): Computes a lower bound on the visibility for which a state becomes unsteerable.

Also, the function [critical_radius](https://github.com/lucporto/LpJm/blob/main/src/critical_radius.jl) is an implementation of the method to study the steerability of two-qubit states proposed in [Nguyen et al. (2019)](https://arxiv.org/abs/1808.09349).

Typical usage can be found in the folder [examples](https://github.com/lucporto/LpJm/tree/main/examples).

The default solver for linear and semidefinite programs is [MOSEK](https://www.mosek.com/), but it can be changed through the kwarg solver.

For more polytopes in dimension three and higher, see [this repository](https://github.com/sebastiendesignolle/ApproximationQuantumStates).
