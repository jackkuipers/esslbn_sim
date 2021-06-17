
local_run <- FALSE
# to run locally we need to set the seed and network size
# this is done externally on the cluster runs
write_data <- FALSE

if (local_run) {
  kk <- 1
  sim_setting <- 1
  seed_number <- 101 # the seed
}

n <- c(20, 80, 140, 200)[kk] # number of nodes
default_alpha <- min(0.4, 20/n) # default value of alpha
N <- 10*n # number of observations
exp_parents <- 2 # expected number of parents
graph_type <- "ERs"
include_truth <- FALSE # don't make artificial search space
MCMC_sample_its <- max(25000, 5*round(n*n*log(n)))

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

if (graph_type == "BA") {
  default_alpha <- min(0.4, 20/n)/2 # default value of alpha
}

# load libraries
library(pcalg)
library(graph)
library(BiDAG)
source("DAGfns.R")

# store values for later dataframe
setup_vec <- c(n, N, exp_parents, graph_type, seed_number)
names(setup_vec) <- c("n", "N", "parents", "type", "seed")
# create a name for the directory to store stuff
f_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
subdir_name <- paste(paste(names(setup_vec[-5]), setup_vec[-5], sep = "_"), collapse = "_")
dir_name <- paste("./sims", subdir_name, sep = "/")


### Generate data
set.seed(seed_number) # set seed
# generate random DAG
if (graph_type == "ERs") {
  trueDAGedges <- as(pcalg::randDAG(n = n, d = 2*exp_parents, 
                                    wFUN = list(runif, min=0.4, max=2)), "matrix")
}
if (graph_type == "ER") {
  trueDAGedges <- as(pcalg::randDAG(n = n, d = 2*exp_parents), "matrix")
}
if (graph_type == "BA") {
  trueDAGedges <- as(pcalg::randDAG(n = n, d = 2*exp_parents, 
                                    method = "barabasi", par1 = 0.4), "matrix")
}
if (graph_type == "iER") {
  trueDAGedges <- as(pcalg::randDAG(n = n, d = 2*exp_parents, 
                                    method = "interEr", par1 = 2, par2 = 0.1), "matrix")
}
trueDAG <- 1*(trueDAGedges != 0)
trueCPDAG <- BiDAG:::dagadj2cpadj(trueDAG)

set.seed(seed_number) # set seed
# generate simulated data
data <- rmvDAG(trueDAGedges, N)

# create directory if none exists
if (!dir.exists("./sims")) { 
  dir.create("./sims")
}
if (!dir.exists(dir_name)) { 
  dir.create(dir_name)
}

if (write_data) {
  if (!file.exists(paste0(dir_name, "/", f_name, "_data.csv"))) {
    write.csv(data, paste0(dir_name, "/", f_name, "_data.csv"), 
                         row.names = FALSE)
  }
}

# to store method name and parameter values
method_vec <- rep(NA, 3)
names(method_vec) <- c("method", "parameter", "value")

### MAP and MCMC search
# fill out missing simulations
if (!file.exists(paste0(dir_name, "/", f_name, "_BiDAG.Rdata"))) {
  
  result_df <- NULL

### MAP search
  am_value <- 0.25

  method_vec[1] <- "MAP"
  method_vec[2] <- "am"
  method_vec[3] <- am_value

  # Start the clock!
  ptm <- proc.time()

  scoreObject <- scoreparameters(scoretype = "bge", data = data, bgepar = list(am = am_value))
  
  set.seed(seed_number) # set seed  
  bestDAGs <- iterativeMCMC(scoreObject, scoreout = TRUE, hardlimit = 16, softlimit = 10,
                            alpha = default_alpha)

  # Stop the clock
  MAPtime <- (proc.time() - ptm)[1]
  
  # Extract CPDAG
  bestCPDAG <-  bestDAGs$CPDAG
  # Compare to generating model
  MAPresult <- compareDAGs(bestCPDAG, trueCPDAG, rnd = 6)

  result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, MAPresult))))

  method_vec[2] <- "iteration"
  for (ii in 1:length(bestDAGs$maxtrace)) {
    method_vec[3] <- ii
    itCPDAG <- BiDAG:::dagadj2cpadj(bestDAGs$maxtrace[[ii]]$DAG)
    itMAPresult <- compareDAGs(itCPDAG, trueCPDAG, rnd = 6)
    result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, itMAPresult))))
  }

## Final order sample

  # Start the clock!
  ptm <- proc.time()

  set.seed(seed_number) # set seed
  orderchain <- orderMCMC(scoreObject, startspace = bestDAGs$endspace, 
                         MAP = FALSE, plus1 = TRUE, chainout = TRUE, 
                         startorder = bestDAGs$maxorder,
                         iterations = MCMC_sample_its,
                         hardlimit = max(colSums(bestDAGs$endspace)))

  # Stop the clock
  MCMCtime <- (proc.time() - ptm)[1]

  method_vec[1] <- "MCMC"
  method_vec[2] <- "p"

  # Compare consensus graphs
  MCMCresult <- samplecomp(orderchain, trueCPDAG, rnd = 6)
  for (ii in 1:nrow(MCMCresult)) {
    method_vec[3] <- as.numeric(MCMCresult[ii, "p"])
    result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, MCMCresult[ii, -9]))))
  }

## Initial order sample

  # Start the clock!
  ptm <- proc.time()

  set.seed(seed_number) # set seed
  iorderchain <- orderMCMC(scoreObject, startspace = bestDAGs$startspace, 
                        MAP = FALSE, plus1 = TRUE, chainout = TRUE, 
                        startorder = bestDAGs$maxorder,
                        iterations = MCMC_sample_its,
                        hardlimit = max(colSums(bestDAGs$startspace)))

  # Stop the clock
  iMCMCtime <- (proc.time() - ptm)[1]

  method_vec[1] <- "iMCMC"
  # Compare consensus graphs
  iMCMCresult <- samplecomp(iorderchain, trueCPDAG, rnd = 6)
  for (ii in 1:nrow(iMCMCresult)) {
    method_vec[3] <- as.numeric(iMCMCresult[ii, "p"])
    result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, iMCMCresult[ii, -9]))))
  }

  if (include_truth) {
  ## Order sample with true DAG included in the search space
    
    joined_space <- bestDAGs$startspace | trueCPDAG
  
  # Start the clock!
  ptm <- proc.time()
  
  set.seed(seed_number) # set seed
  morderchain <- orderMCMC(scoreObject, startspace = joined_space, 
                           MAP = FALSE, plus1 = TRUE, chainout = TRUE, 
                           startorder = bestDAGs$maxorder,
                           iterations = MCMC_sample_its,
                           hardlimit = max(colSums(joined_space)))
  
  # Stop the clock
  mMCMCtime <- (proc.time() - ptm)[1]
  
  method_vec[1] <- "mMCMC"
  # Compare consensus graphs
  mMCMCresult <- samplecomp(morderchain, trueCPDAG, rnd = 6)
  for (ii in 1:nrow(mMCMCresult)) {
    method_vec[3] <- as.numeric(mMCMCresult[ii, "p"])
    result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, mMCMCresult[ii, -9]))))
  }
  
  }
  
  time_df <- data.frame(n = n, N = N, parents = exp_parents, 
                        seed = seed_number,
                      method = c("MAP", "MCMC", "iMCMC"), 
                      time = c(MAPtime, MCMCtime, iMCMCtime))

  save(result_df, time_df, file = paste0(dir_name, "/", f_name, "_BiDAG.Rdata"))
}

## GES
# fill out missing simulations
if (!file.exists(paste0(dir_name, "/", f_name, "_GES.Rdata"))) {
  
  result_df <- NULL  
  
  # Start the clock!
  ptm <- proc.time()
  
  method_vec[1:2] <- c("GES", "lambda")
  ges_lambdas <- c(1,2,5,7,9,11,15,25) 

  for (lambda_scale in ges_lambdas) {
    
    method_vec[3] <- lambda_scale
    ges_score <- new("GaussL0penObsScore", data, 
                      lambda = lambda_scale*log(nrow(data)))
    
    set.seed(seed_number) # set seed
    ges_fit <- ges(ges_score)
    
    # Extract CPDAG
    gesCPDAG <- 1*as(ges_fit$essgraph, "matrix")
    # Compare to generating model
    GESresult <- compareDAGs(gesCPDAG, trueCPDAG, rnd = 6)
    
    result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, GESresult))))
  }
  
  # Stop the clock - this is over all penalisations
  GEStime <- (proc.time() - ptm)[1]
  
  time_df <- data.frame(n = n, N = N, parents = exp_parents, 
                        seed = seed_number,
                        method = c("GES"), 
                        time = c(GEStime))
  
  save(result_df, time_df, file = paste0(dir_name, "/", f_name, "_GES.Rdata"))
}

## PC
# fill out missing simulations
if (!file.exists(paste0(dir_name, "/", f_name, "_PC.Rdata"))) {
  
  result_df <- NULL 

  # correlation matrix  
  cor_mat <- cor(data)
  
  # Start the clock!
  ptm <- proc.time()
  
  method_vec[1:2] <- c("PC", "alpha")
  pc_alphas <- c(0.01, 0.05, 0.1, 0.2, 0.35, 0.45)
  
  if (sim_setting == 5 && n > 100){ # truncate as it takes too long otherwise
    pc_alphas <- pc_alphas[1:3]
  }

  for (alpha in pc_alphas) {
    
    method_vec[3] <- alpha
    
    set.seed(seed_number) # set seed
    ## estimate CPDAG
    pc_fit <- pc(suffStat = list(C = cor_mat, n = nrow(data)),
                 indepTest = gaussCItest, ## indep.test: partial correlations
                 alpha = alpha, labels = paste0("V", 1:ncol(data)))
    # Extract CPDAG
    pcCPDAG <- 1*as(pc_fit@graph, "matrix")
    # Compare to generating model
    PCresult <- compareDAGs(pcCPDAG, trueCPDAG, rnd = 6)
    
    result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, PCresult))))
  }
  
  # Stop the clock - over all alphas
  PCtime <- (proc.time() - ptm)[1]
  
  time_df <- data.frame(n = n, N = N, parents = exp_parents, 
                        seed = seed_number,
                        method = c("PC"), 
                        time = c(PCtime))
  
  save(result_df, time_df, file = paste0(dir_name, "/", f_name, "_PC.Rdata"))
}

## PC time for initialisation
# fill out missing simulations
if (!file.exists(paste0(dir_name, "/", f_name, "_PCt.Rdata"))) {
  
  # correlation matrix  
  cor_mat <- cor(data)
  
  # Start the clock!
  ptm <- proc.time()
  
  method_vec[1:2] <- c("PC", "alpha")
  alpha <- default_alpha 
  method_vec[3] <- alpha
    
  set.seed(seed_number) # set seed
  ## estimate CPDAG
  pc_fit <- pc(suffStat = list(C = cor_mat, n = nrow(data)),
                 indepTest = gaussCItest, ## indep.test: partial correlations
                 alpha = alpha, labels = paste0("V", 1:ncol(data)))
  # Extract CPDAG
  pcCPDAG <- 1*as(pc_fit@graph, "matrix")

  # Stop the clock - over all alphas
  PCttime <- (proc.time() - ptm)[1]
  
  time_df <- data.frame(n = n, N = N, parents = exp_parents, 
                        seed = seed_number,
                        method = c("PCt"), 
                        time = c(PCttime))
  
  save(time_df, file = paste0(dir_name, "/", f_name, "_PCt.Rdata"))
}

