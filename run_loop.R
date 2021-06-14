for (kk in 1:4) {
  for (jj in 1:6) {
    sim_setting <- jj
    for (seedy in 101:200) {
      seed_number <- seedy
      print(paste0("setting", jj, "size", kk, "seed", seedy))
      source("./esslbn_sim.R")
    }
  }
}
