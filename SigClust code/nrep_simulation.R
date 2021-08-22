rm(list=ls())
setwd('C:/Users/Steve/OneDrive - Imperial College London/MSc/Summer Project/R Code/git_simulations_and_experiments')
# install.packages("sigclust")
library(sigclust) # for sigclust()
library(MASS) # for mvrnorm() 
library(ggplot2) # for ggplot() 
library(tictoc) # for tic() and toc()
library(reshape2) # for melt()
library(foreach) # for foreach()
library(parallel) # for detectCores()
library(doParallel) # for registerDoParallel()

numCores = detectCores()
registerDoParallel(numCores)

my_theme = theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12))

### CASE OF ONE CLUSTER
###-------------------------------------------------------------------------- 
n = 40
d = 80
mu = rep(0, d)
D = diag(d)
v = c(1, 64, 81, 100) # different values for first eigenvalue
D[1, 1] = v[4] # v = 1
D[1, ]

# using parallel processing is much faster - had checked with tic() toc()
pvalue1 = foreach (b=1:100) %dopar%{
  library(MASS) # for mvrnorm() 
  library(sigclust) # for sigclust()
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=56 d-dimensional observations
  nsim = 1000
  nrep = 1
  icovest = 3
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}
pval1 = sapply(1:100, function(i) pvalue1[[i]]@pval)

pvalue2 = foreach (b=1:100) %dopar%{
  library(MASS) # for mvrnorm() 
  library(sigclust) # for sigclust()
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=56 d-dimensional observations
  nsim = 1000
  nrep = 2
  icovest = 3
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}
pval2 = sapply(1:100, function(i) pvalue2[[i]]@pval)

pvalue3 = foreach (b=1:100) %dopar%{
  library(MASS) # for mvrnorm() 
  library(sigclust) # for sigclust()
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=56 d-dimensional observations
  nsim = 1000
  nrep = 3
  icovest = 3
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}
pval3 = sapply(1:100, function(i) pvalue3[[i]]@pval)

save(pval1, pval2, pval3, file="pvals_nrep_sim.RData")

load("pvals_nrep_sim.RData")

pvals = as.data.frame(cbind(pval1, pval2, pval3))
names(pvals) = c(1, 2, 3)
pvals = melt(pvals, measure.vars=c("1", "2", "3"), variable.name="nrep", 
             value.name="pvalue")
dev.off()

ggplot() + stat_ecdf(dat=pvals, aes(x=pvalue, linetype=nrep, colour=nrep)) + 
  xlim(0, 1) + geom_vline(xintercept=0.05) +
  labs(x="p-value", y="Empirical Distribution") + 
  ggtitle("Empirical cdf's of SigClust p values for different values of nrep") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate(geom = "text", label = bquote(alpha == .(0.05)), x = 0.1, y = 0.9, 
           size=3.2) + 
  scale_colour_brewer(palette="Set1", direction=-1) +
  my_theme
# therefore nrep does not seem to be having any effect on the obtained pvalues
# and subsequent ecdf - so we are safe to use nrep=1 to optimise speed