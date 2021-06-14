# esslbn_sim
Simulation scripts to benchmark sampling and structure learning of Bayesian networks

To accompany https://arxiv.org/abs/1803.07859

The file esslbn_sim.R runs simulations for a single setting, seed and network size. This may be run locally, or called from an outside loop (like run_loop.R) for example on a cluster.

The file collateResults.R collates output files from the different runs into the sims_collated directory.
