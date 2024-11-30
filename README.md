# Group-driven Voter Model (GVM) on Hypergraphs 

This repository contains the codes to reproduce Monte Carlo (MC) simulations associated with the paper:
"Competition between group interactions and nonlinearity in voter dynamics on hypergraphs" by Jihye Kim, Deok-Sun Lee, Byungjoon Min, Mason A. Porter, Maxi San Miguel, and K.-I. Goh (2024). 

## Abstract

Social dynamics are often driven by both pairwise (i.e., dyadic) relationships and higher-order (i.e., polyadic) group relationships, which one can describe using hypergraphs. To gain insight into the impact of polyadic relationships on dynamical processes on networks, we formulate and study a polyadic voter process, which we call the group-driven voter model (GVM), in which we incorporate the effect of group interactions by nonlinear interactions that are subject to a group (i.e., hyperedge) constraint. By examining the competition between nonlinearity and group sizes, we show that the GVM achieves consensus faster than standard voter-model dynamics, with an optimum minimizing exit time $\tau$. We substantiate this finding by using mean-field theory on annealed uniform hypergraphs with $N$ nodes, for which $\tau$ scales as $\mathcal{A} \ln N$, where the prefactor $\mathcal{A}$ depends both on the nonlinearity and on group-constraint factors. Our results reveal how competition between group interactions and nonlinearity shapes GVM dynamics. We thereby highlight the importance of such competing effects in complex systems with polyadic interactions.

## MC simulation codes

#### 1. Fig2_geometric_Ps.cpp
This code is used for MC simulations for simplicial GVM on annealed hypergraphs with $N$ nodes and geometric hyperedge-size distribution. 
One can reproduce Fig. 2(a) of the main manuscript by using this code.

#### 2. Fig2_powerlaw_Ps.cpp
This code is used for MC simulations for simplicial GVM on annealed hypergraphs with $N$ nodes and power law hyperedge-size distribution.
One can reproduce Fig. 2(b) of the main manuscript by using this code.

#### 3. Fig3_complete_networks.cpp 
This code is used for MC simulations for our GVM on complete networks with $N$ nodes.
One can reproduce the blue data points in the inset of Fig. 3 of the main manuscript by using this code.

#### 4. Fig3_with_small_s.cpp
This code is used for MC simulations for our GVM on annealed s-uniform hypergraphs with $N$ nodes.
One can reproduce the data points (for the case of $s=3$) of Fig. 3 of the main manuscript by using this code.

#### 5. Fig4_complete_networks.cpp 
This code is used for MC simulations for our GVM on complete networks with $N$ nodes.
One can reproduce the blue data points of Fig. 4 of the main manuscript by using this code.

#### 6. Fig4_with_small_s.cpp
This code is used for MC simulations for our GVM on annealed s-uniform hypergraphs with $N$ nodes.
One can reproduce the red and green data points (for the cases of $s=3$ and $s=7$, respectively) of Fig. 4 of the main manuscript by using this code.

#### 7. Fig5a_complete_networks.cpp 
This code is used for MC simulations for our GVM on complete networks with $N$ nodes.
One can reproduce the black data points of Fig. 5(a) of the main manuscript by using this code.

#### 8. Fig5a_with_small_s.cpp
This code is used for MC simulations for our GVM on annealed s-uniform hypergraphs with $N$ nodes.
One can reproduce the data points (for the cases of $s \ll N$) of Fig. 5(a) of the main manuscript by using this code.

#### 9. FigS2_initial_density_dependence.cpp
This code is used for MC simulations for our GVM on annealed s-uniform hypergraphs with $N$ nodes.
One can reproduce the data points of Fig. S2 of the Supplementary Material (SM) by using this code.

#### 9. FigS3_geometric_Ps.cpp
This code is used for MC simulations for our GVM on annealed hypergraphs with $N$ nodes and geometric hyperedge-size distribution.
One can reproduce the data points of Fig. S3 of the Supplementary Material (SM) by using this code.

#### 9. FigS4_powerlaw_Pk.cpp
This code is used for MC simulations for our GVM on annealed hypergraphs with $N$ nodes and powerlaw degree distribution.
One can reproduce the data points of Fig. S4 of the Supplementary Material (SM) by using this code.

#### 10. FigS5_prohibit_duplicate_selection.cpp
This code is used for MC simulations for GVM without duplicate selections on annealed s-uniform hypergraphs with $N$ nodes.
One can reproduce the red and green data points (for the cases of $s=3$ and $s=7$, respectively) of Fig. S3 of the SM by using this code.

#### 11. FigS5_prohibit_duplicate_selection_on_complete_networks.cpp
This code is used for MC simulations for GVM without duplicate selections on complete networks with $N$ nodes.
One can reproduce the blue data points of Fig. S3 of the SM by using this code.

#### 12. FigS6_hyperedge_flipping.cpp
This code is used for MC simulations for GVM with hyperedge-update dynamics on annealed s-uniform hypergraphs with $N$ nodes.
One can reproduce the data points of Fig. S4 of the SM by using this code.
