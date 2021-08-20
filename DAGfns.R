
### This function gives edges weights between the bounds
# with both positive and negative signs

wFUN <- function(m, lb, ub){ # function for edge weights
  runif(m, lb, ub)*sample(c(-1, 1), m, replace = TRUE)
}


### This function generates Gaussian data from a DAG
# following the topological order

rmvDAG <- function(trueDAGedges, N, standardise = TRUE) {
  trueDAG <- 1*(trueDAGedges != 0) # the edge presence in the DAG
  n <- ncol(trueDAG) # number of variables
  data <- matrix(0, nrow = N, ncol = n) # to store the simulated data
  top_order <- rev(BiDAG:::DAGtopartition(n, trueDAG)$permy) # go down order
  for (jj in top_order) {
    parents <- which(trueDAG[, jj] == 1) # find parents
    lp <- length(parents) # number of parents
    if (lp == 0) { # no parents
      data[, jj] <- 0
    } else if (lp == 1) { # one parent
      data[, jj] <- data[, parents]*trueDAGedges[parents, jj]
    } else { # more than one parent
      data[, jj] <- colSums(t(data[, parents])*trueDAGedges[parents, jj])
    }
    # add random noise
    data[, jj] <- data[, jj] + rnorm(N)
  }
  if(standardise) { # whether to standardise
    scale(data)
  } else {
    data
  }
}

samplecomp <- function (MCMCchain, truedag, p = c(0.99, 0.95, 0.9, 0.8, 0.7, 
                                                  0.6, 0.5, 0.4, 0.3, 0.2), pdag = TRUE, burnin = 0.2, trans = TRUE, rnd = 2) 
{
  if (is.matrix(truedag)) 
    truedag <- m2graph(truedag)
  MCMCmatlist <- MCMCchain$traceadd$incidence
  n <- nrow(MCMCmatlist[[1]])
  truecp <- pcalg::dag2cpdag(truedag)
  if (MCMCchain$info$DBN) {
    pdag <- FALSE
    if (trans == TRUE) {
      cat("comparison is performed for transition structures \n")
      trueadj <- DBNcut(graph2m(truedag), MCMCchain$info$nsmall, 
                        MCMCchain$info$bgn)
      truedag <- m2graph(trueadj)
      truecpadj <- graph2m(truecp)
      truecpadj <- DBNcut(truecpadj, MCMCchain$info$nsmall, 
                          MCMCchain$info$bgn)
      truecp <- m2graph(truecpadj)
    }
    else {
      cat("comparison is performed for initial structures \n")
      trueadj <- DBNinit(graph2m(truedag), MCMCchain$info$nsmall, 
                         MCMCchain$info$bgn)
      truedag <- m2graph(trueadj)
      truecpadj <- DBNinit(graph2m(truecp), MCMCchain$info$nsmall, 
                           MCMCchain$info$bgn)
      truecp <- m2graph(truecpadj)
      n <- MCMCchain$info$nsmall + MCMCchain$info$bgn
    }
  }
  endstep <- length(MCMCmatlist)
  startstep <- max(as.integer(burnin * endstep), 1)
  if (pdag) {
    dags <- lapply(MCMCmatlist[startstep:endstep], BiDAG:::dagadj2cpadj)
  }
  else {
    dags <- MCMCmatlist[startstep:endstep]
  }
  if (MCMCchain$info$DBN) {
    if (trans) {
      dags <- lapply(dags, DBNcut, n.dynamic = MCMCchain$info$nsmall, 
                     n.static = MCMCchain$info$bgn)
    }
    else {
      dags <- lapply(dags, DBNinit, n.dynamic = MCMCchain$info$nsmall, 
                     n.static = MCMCchain$info$bgn)
    }
  }
  postprobmat <- Reduce("+", dags)/(endstep - startstep + 1)
  if (length(p) == 1) {
    mlist <- matrix(0, nrow = n, ncol = n)
    mlist[which(postprobmat > p)] <- 1
    res <- compareDAGs(mlist, truedag)
    res <- c(res, p)
    names(res)[length(res)] <- "p"
  }
  else {
    mlist <- list()
    i <- 1
    for (py in 1:length(p)) {
      mlist[[i]] <- matrix(0, nrow = n, ncol = n)
      mlist[[i]][which(postprobmat > p[py])] <- 1
      i <- i + 1
    }
    res <- lapply(mlist, compareDAGs, truedag, rnd = rnd)
    res <- Reduce(rbind, res)
    res <- cbind(res, p)
    rownames(res) <- c(1:nrow(res))
  }
  return(res)
}
