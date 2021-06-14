
### range of seeds
seed_numbers <- 100 + 1:100 # the seeds

### Settings
for (n in c(20, 80, 140, 200)) { # number of nodes
  for (sim_setting in 1:6){ # different simulation settings
    N <- 10*n # number of observations
    exp_parents <- 2 # expected number of parents
    graph_type <- "ERs"
    include_truth <- FALSE
    
    if (sim_setting == 2){ # smaller sample size
      N <- 2*n # number of observations
    }
    if (sim_setting == 3){ # ER graphs
      graph_type <- "ER"
    }
    if (sim_setting == 4){ # BA graphs
      graph_type <- "BA"
    }
    if (sim_setting == 5){ # BA graphs, fewer parents
      graph_type <- "BA"
      exp_parents <- 1 # expected number of parents
    }
    if (sim_setting == 6){ # iER graphs, fewer parents
      graph_type <- "iER"
      exp_parents <- 3 # expected number of parents
    }
    

# store values for later dataframe
setup_vec <- c(n, N, exp_parents, graph_type)
names(setup_vec) <- c("n", "N", "parents", "type")
# create a name for the directory to store stuff
subdir_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
dir_name <- paste("./sims", subdir_name, sep = "/")
dir_name2 <- paste("./sims_collated", subdir_name, sep = "/")

method_files <- c("BiDAG", "GES", "PC", "PCt")

#if (!file.exists(paste0(dir_name, ".Rdata"))) {

results_df <- NULL
times_df <- NULL

for (seed_number in seed_numbers) {
  seed_part <- paste("", "seed", seed_number, sep = "_")
  for (ii in 1:length(method_files)) {
    if (file.exists(paste0(dir_name, "/", subdir_name, seed_part, "_", method_files[ii], ".Rdata"))) {
      load(paste0(dir_name, "/", subdir_name, seed_part, "_", method_files[ii], ".Rdata"))
      results_df <- rbind(results_df, result_df)
      times_df <- rbind(times_df, time_df)
    }
  }
}

save(results_df, times_df, file = paste0(dir_name2, ".Rdata"))

  }
}
