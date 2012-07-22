Ranked gene tree topologies probability computation
=================
A polynomial-time algorithm has been described for computing probabilities of ranked gene tree topologies given species trees. Once, implemented, ranked gene tree probabilities could be used to infer species trees, although inferring species trees is beyond the scope of the project. The idea is to consider ranked gene tree topologies, where we distinguish the relative order of times of nodes on gene trees, but not the real-valued branch lengths.


This project will provide an implementation of polynomial time algorithm to calculate the probability of a ranked gene tree topology for a given species tree.

TODO:
indexed speciation array should contain distance information to get rid of future get_dist_from_root calls:
node2int --> node2float

Indexing:

s[] (speciation): same as in paper
u[] (coalescence): 0-based (index as u[i-1])
m[] (intervals in tau[i]): 0-based (index as m[i-1])
gene_lineages (L[]): 0-based (index as gene_lineages[i-1])

tau_idx indiсes are calculated from s[] (compatible)

y[][] (beaded tree index): 0-based (index as y[i-1][j-1]) ?
g[] array:
k[][][] array:

LCA:
 *  equivalent node ids are fetched with taxa, and use s[] indices


indexed species array is indexed with s[] indices

You will need to install libgmp to use this code (apt-get install libgmp-dev)
