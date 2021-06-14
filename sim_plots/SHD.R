library(gridExtra)
library(grid)
source(file="./gridarrange.R")

library(tidyverse)

medSHD_df <- NULL
all_results_df <- NULL

for(sim_setting in 1:2){
  
  pl <- vector("list", 4)
  
  for (kk in 1:4) {

n <- c(20, 80, 140, 200)[kk]
N <- 10*n # number of observations
exp_parents <- 2 # expected number of parents
graph_type <- "ERs"

if (sim_setting == 2){ # smaller sample size
  N <- 2*n # number of observations
}

load(paste0("../sims_collated/n_", n, "_N_", N, "_parents_",
            exp_parents, "_type_", graph_type, ".Rdata"))

results_df %>% filter(method != "mMCMC", parameter !="iteration") %>%
  mutate(SHDnum = as.numeric(as.character(SHD)),
         nnum = as.numeric(as.character(n)),
         Nnum = as.numeric(as.character(N)),
         val = as.numeric(as.character(value)),
         SHDn = SHDnum/nnum) -> results_df_processed

results_df_processed %>% group_by(n, N, method, parameter, val) %>% 
  summarise(mSHDn = median(SHDn)) -> summar_df

medSHD_df <- rbind(medSHD_df, summar_df)
all_results_df <- rbind(all_results_df, results_df_processed)

  }
}

medSHD_df %>% group_by(method, parameter, val) %>% 
  summarise(meanSHDn = mean(mSHDn)) -> summy_df

summy_df %>% group_by(method) %>% summarise(val = val[which.min(meanSHDn)])

# PC has same value at 0.01 and 0.05 so we go with default at 0.05
# GES has best value at 2
# MCMC best around 0.5 or 0.6 we take 0.6

all_results_df %>% filter(val %in% c(0.25, 0.6, 2, 0.05)) %>% 
  select(-n, -N) %>% mutate(algorithm = method) -> SHD_df

SHD_df %>% mutate(ss = Nnum/nnum) %>% group_by(ss) %>% summarise(maxSHDn = max(SHDn))

SHD_df$algorithm = with(SHD_df, factor(algorithm, 
                        levels=c("iMCMC", "MAP", "MCMC", "PC", "GES")))

for(sim_setting in 1:2){
  
  pl <- vector("list", 4)
  
  for (kk in 1:4) {
    
    n <- c(20, 80, 140, 200)[kk]
    N <- 10*n # number of observations
    if (sim_setting == 2){ # smaller sample size
      N <- 2*n # number of observations
    }
    graph_title <- paste0("n=", n, ", sample size ",N/n,"n") 

pl[[kk]] <- ggplot(SHD_df %>% filter(nnum == n, Nnum == N),
             aes(x = algorithm, y = SHDn, col = algorithm)) +
  geom_boxplot(aes(fill = algorithm)) + theme(axis.title.x = element_blank()) +
  scale_color_manual(labels=c('MCMC\ninitial','MCMC\nMAP','MCMC\nfinal','PC','GES'),
                    values = c("pink", "red", "purple",  "#377eb8", "#4daf4a"))+
  scale_fill_manual(labels=c('MCMC\ninitial','MCMC\nMAP','MCMC\nfinal','PC','GES'),
                     values =c("#f4cae4", "#fb8072", "#bebada",  "#a6cee3", "#b2df8a") ) +
  scale_x_discrete(labels = c('MCMC\ninitial','MCMC\nMAP',
                              'MCMC\nfinal','PC','GES')) +
  ylim(c(0,4.5))+ggtitle(graph_title)
}

  pdf(file = paste0("./SHD_N_", N/n, "n_parents_",
                    exp_parents, "_type_", graph_type, ".pdf"), width=6.5,height=4.2, onefile=FALSE)
  
  grid_arrange_shared_legend(pl[[1]],  pl[[2]],  pl[[3]],  pl[[4]], ncol = 2, nrow = 2,position="right")
  
  dev.off()
  
}
