# esslbn_sim
Simulation scripts to benchmark sampling and structure learning of Bayesian networks

To accompany https://arxiv.org/abs/1803.07859

The file esslbn_sim.R runs simulations for a single setting, seed and network size. This may be run locally, or called from an outside loop (like run_loop.R) for example on a cluster.

The file collateResults.R collates output files from the different runs into the sims_collated directory.

# session info

The files have been tested with the following (summer 2021) package versions: 
* R version 4.0.2
* BiDAG version 2.0.3
* pcalg version 2.7-3
* graph version 1.66.0
