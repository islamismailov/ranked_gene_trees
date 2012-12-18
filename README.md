Ranked gene tree topologies probability computation
=================
This is a polynomial-time algorithm implementation for computing probabilities of ranked gene tree topologies given species trees. Ranked gene tree probabilities could be used to infer species trees, although inferring species trees is beyond the scope of the project. The idea is to consider ranked gene tree topologies, where we distinguish the relative order of times of nodes on gene trees, but not the real-valued branch lengths.

This project provides an efficient implementation of polynomial time algorithm to calculate the probability of a ranked gene tree topology for a given species tree.

You will need to install GMP and MPFR to use this code (apt-get install libgmp-dev libmpfr-dev)

