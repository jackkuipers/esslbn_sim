library(gridExtra)
library(grid)
source(file="./gridarrange.R")

library(tidyverse)

for(sim_setting in 1:6){

pl <- vector("list", 4)

for (kk in 1:4) {

n <- c(20, 80, 140, 200)[kk]
N <- 10*n # number of observations
exp_parents <- 2 # expected number of parents
graph_type <- "ERs"

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
  include_truth <- TRUE # make artificial search space
}

if(graph_type == "ERs"){
  graph_title <- paste0("n=", n, ", sample size ",N/n,"n") 
} else {
  graph_title <- bquote(plain("n=")*.(n)*plain(", N=")*.(N/n)*plain("n, ")*nu*plain("=")*.(exp_parents))
}


load(paste0("../sims_collated/n_", n, "_N_", N, "_parents_",
            exp_parents, "_type_", graph_type, ".Rdata"))

results_df %>% mutate(TPRnum = as.numeric(as.character(TPR)),
                      FPRnnum = as.numeric(as.character(FPRn)),
                      val = as.numeric(as.character(value))) %>% 
  group_by(method, parameter, val) %>% filter(method != "mMCMC") %>%
  summarise(TPR = median(TPRnum),
            q1TPR = quantile(TPRnum, 0.25),
            q3TPR = quantile(TPRnum, 0.75),
            FPRp = median(FPRnnum), n = n()) %>%
  mutate(algorithm = method) -> summar_df 

summar_df$algorithm = with(summar_df, factor(algorithm, 
                    levels=c("GES", "PC", "iMCMC", "MCMC", "MAP")))

width.er<-0.01

ggplot(data = summar_df %>% filter(parameter != "iteration") %>% 
         filter(algorithm != "MAP"), 
  aes(x = FPRp, y = TPR, group = val, col = algorithm)) + 
  geom_errorbar(aes(ymin = q1TPR, ymax = q3TPR, col = algorithm), width = width.er) + 
  geom_path(aes(group = algorithm, col = algorithm)) +
  geom_point(aes(group = algorithm, col = algorithm, shape = algorithm), size = 2) +
  coord_cartesian(xlim = c(-1.25*width.er, 0.4),
    ylim = c(0, 1 + 5*width.er), expand = FALSE) + # zoom plot rather than truncate
  #xlim(c(-width.er, 0.4 + 0.1*(sim_setting>2) + width.er)) + ylim(c(0, 1)) +
  scale_shape_manual(labels=c("GES","PC","MCMC\ninitial","MCMC\nfinal"), values =c(8,17,16,15))+
  scale_colour_manual(labels=c("GES","PC","MCMC\ninitial","MCMC\nfinal"), values = c("#4daf4a","#377eb8","pink","purple")) +
  ggtitle(graph_title) -> pl[[kk]]

}

#sample size 10n

#pdf(file="ROC10n.pdf",width=6.5,height=4.2, onefile=FALSE)

pdf(file = paste0("./ROC_N_", N/n, "n_parents_",
            exp_parents, "_type_", graph_type, ".pdf"), width=6.5,height=4.2, onefile=FALSE)

grid_arrange_shared_legend(pl[[1]],  pl[[2]],  pl[[3]],  pl[[4]], ncol = 2, nrow = 2,position="right")

dev.off()

}
