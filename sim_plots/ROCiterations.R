library(gridExtra)
library(grid)
source(file="./gridarrange.R")

library(tidyverse)

pl <- vector("list", 4)

for (kk in 1:4) {

n <- c(20, 80, 140, 200)[kk]
N <- 10*n # number of observations
exp_parents <- 2 # expected number of parents
graph_type <- "ERs"

graph_title <- paste0("n=", n, ", sample size ",N/n,"n") 

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
  filter(n > 90) %>%
  mutate(algorithm = method) %>% 
  mutate(val2 = val + 100*(parameter=="am"))-> summar_df 

summar_df$algorithm = with(summar_df, factor(algorithm, 
                    levels=c("GES", "PC", "iMCMC", "MAP", "MCMC")))

summar_df2 <- rbind(summar_df %>% filter(parameter %in% c("iteration", "p", "am")) %>%
  filter(parameter != "am"),
  summar_df %>% filter(parameter %in% c("iteration", "p", "am")) %>%
    filter(parameter == "am")
)

width.er<-0.01

ggplot(data = summar_df2, 
  aes(x = FPRp, y = TPR, group = as.factor(val2), col = algorithm)) + 
  #geom_errorbar(aes(ymin = q1TPR, ymax = q3TPR, col = algorithm), width = width.er) + 
  geom_path(aes(group = algorithm, col = algorithm, linetype = algorithm)) +
  geom_point(aes(group = algorithm, col = algorithm, shape = algorithm), size = 2) +
  coord_cartesian(xlim = c(-1*width.er, 0.18 + width.er),
    ylim = c(0.6 - 2*width.er, 1 + 2*width.er), expand = FALSE) + # zoom plot rather than truncate
  #xlim(c(-width.er, 0.4 + 0.1*(sim_setting>2) + width.er)) + ylim(c(0, 1)) +
  scale_linetype_manual(labels=c("MCMC\ninitial","MAP\niterations","MCMC\nfinal"),values = c(1, 3,1))+
  scale_shape_manual(labels=c("MCMC\ninitial","MAP\niterations","MCMC\nfinal"),values = c(16, 3,15))+
  scale_colour_manual(labels=c("MCMC\ninitial","MAP\niterations","MCMC\nfinal"),values = c("pink", "red","purple")) +
  ggtitle(graph_title) -> pl[[kk]]

}

pdf(file = paste0("./ROC_iterations.pdf"), width=6.5,height=4.2, onefile=FALSE)

grid_arrange_shared_legend(pl[[1]],  pl[[2]],  pl[[3]],  pl[[4]], ncol = 2, nrow = 2,position="right")

dev.off()

