library(gridExtra)
library(grid)
source(file="./gridarrange.R")

library(tidyverse)
options(encoding = "UTF-8")

pl <- vector("list", 4)

shifty <- 2

for(sim_setting in 3:6){

  mtime_df <- NULL
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
      graph_title <- expression("Erdös-Rényi, "*nu*"=2")
    }
    if (sim_setting == 4){ # BA graphs
      graph_type <- "BA"
      graph_title <- expression("Bárabasi-Albert, "*nu*"=2")
    }
    if (sim_setting == 5){ # BA graphs, fewer parents
      graph_type <- "BA"
      exp_parents <- 1 # expected number of parents
      graph_title <- expression("Bárabasi-Albert, "*nu*"=1")
    }
    if (sim_setting == 6){ # iER graphs, fewer parents
      graph_type <- "iER"
      exp_parents <- 3 # expected number of parents
      graph_title <- expression("Erdös-Rényi islands, "*nu*"=3")
    }
    
    load(paste0("../sims_collated/n_", n, "_N_", N, "_parents_",
                exp_parents, "_type_", graph_type, ".Rdata"))

    times_df %>% filter(method %in% c("GES", "PCt", "iMCMC", "MAP", "MCMC")) %>%
      group_by(n, N, parents, method) %>% # GEST time is for all parameters in the ROC curve!
      summarise(mtime_temp = mean(time)*10^shifty) %>% 
      mutate(mtime = mtime_temp/(1+7*(method=="GES"))) -> time_df

    mtime_df <- rbind(mtime_df, time_df)
  }

mtime_df$algorithm <- with(mtime_df,factor(method,levels=c("GES", "PCt", "iMCMC", "MAP","MCMC")))

  
pl[[sim_setting - 2]] <- ggplot(data=mtime_df, aes(x = as.factor(n), y = mtime, group=algorithm)) +   
  labs(y = "runtime in seconds", x = "n") +
  labs(fill = "") +
  scale_fill_manual(labels=c("\nGES\n","\nPC\n","\nMCMC\ninitial\nsample\n","\nMCMC\niterative\nsearch\n","\nMCMC\nsample\n"),
                    values = alpha(c("#4daf4a","#377eb8","pink", "red","purple"),0.6)) +
  scale_y_continuous(trans ="log10", breaks = 10^(0:3+shifty), labels = 10^(0:3))+
  coord_cartesian(ylim = c(0.03, 3000)*10^shifty) + 
  geom_bar(aes(fill = algorithm), position="dodge",  stat = 'identity') +
  ggtitle(graph_title)# + theme(legend.key.size = unit(0.8, 'cm'))
}

pdf(file="Runtime.pdf",width=6.5,height=4.2,  onefile=FALSE)

grid_arrange_shared_legend(pl[[1]], pl[[3]], pl[[2]], pl[[4]], ncol = 2, nrow = 2,position="right")

dev.off()
