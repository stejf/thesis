rm(list=ls())
library(sigclust) # for sigclust()
library(MASS) # for mvrnorm() 
library(tictoc) # for tic() and toc()
library(foreach) # for foreach()
library(parallel) # for detectCores()
library(doParallel) # for registerDoParallel()

numCores = detectCores()
registerDoParallel(numCores)

### CASE OF MIXTURE OF TWO GAUSSIAN DISTRIBUTIONS WITH SIGNAL IN ONE DIRECTION 
###-------------------------------------------------------------------------- 
n = 40
d = 300
mu = rep(0, d)
a = c(0, 10, 20, 25)
D = diag(d)
v = 35
D[1, 1] = v
D[1, ]

pvalue1 = foreach (b=1:100, .packages=c('MASS', 'sigclust')) %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=40 d-dimensional observations
  dat[(n/2+1):(n), 1] = dat[(n/2+1):(n), 1] + a[1]
  nsim = 1000
  nrep = 1
  icovest = 3
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}
###-------------------------------------------------------------------------- 

pvalue2 = foreach (b=1:100, .packages=c('MASS', 'sigclust')) %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=40 d-dimensional observations
  dat[(n/2+1):(n), 1] = dat[(n/2+1):(n), 1] + a[2]
  nsim = 1000
  nrep = 1
  icovest = 3
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}

###-------------------------------------------------------------------------- 
pvalue3 = foreach (b=1:100, .packages=c('MASS', 'sigclust')) %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=40 d-dimensional observations
  dat[(n/2+1):(n), 1] = dat[(n/2+1):(n), 1] + a[3]
  nsim = 1000
  nrep = 1
  icovest = 3
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}

###-------------------------------------------------------------------------- 
pvalue4 = foreach (b=1:100, .packages=c('MASS', 'sigclust')) %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=40 d-dimensional observations
  dat[(n/2+1):(n), 1] = dat[(n/2+1):(n), 1] + a[4]
  nsim = 1000
  nrep = 1
  icovest = 3
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}
save(pvalue1, pvalue2, pvalue3, pvalue4, file="pvalues.RData")
###-------------------------------------------------------------------------- 
#Plots
rm(list=ls())
setwd('C:/Users/Steve/OneDrive - Imperial College London/MSc/Summer Project/R Code/Huang d_300/3_2hard')
load("pvalues.RData")
library(ggplot2) # for ggplot() 
library(reshape2) # for melt()

my_theme = theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12))

a = c(0, 10, 20, 25)

pval1 = sapply(1:100, function(i) pvalue1[[i]]@pval)
pval2 = sapply(1:100, function(i) pvalue2[[i]]@pval)
pval3 = sapply(1:100, function(i) pvalue3[[i]]@pval)
pval4 = sapply(1:100, function(i) pvalue4[[i]]@pval)

# mean p-values
mean(pval1)
mean(pval2)
mean(pval3)
mean(pval4)

# number of p-values less than 0.05
sum(pval1 < 0.05)
sum(pval2 < 0.05)
sum(pval3 < 0.05)
sum(pval4 < 0.05)

# number of p-values less than 0.1
sum(pval1 < 0.1)
sum(pval2 < 0.1)
sum(pval3 < 0.1)
sum(pval4 < 0.1)

pvals = as.data.frame(cbind(pval1, pval2, pval3, pval4))
names(pvals) = a
pvals = melt(pvals, measure.vars=as.character(a), variable.name="a", 
             value.name="pvalue")

p = ggplot() + stat_ecdf(dat=pvals, aes(x=pvalue, linetype=a, colour=a)) + 
  xlim(0, 1) + geom_vline(xintercept=0.05, color="gray45") +
  labs(x="p-value", y="Empirical Distribution") + 
  ggtitle("Mixture of Gaussians - Signal in 1 Direction (hard)") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate(geom = "text", label = bquote(alpha == .(0.05)), x = 0.05, y = 1.1, 
           size=5) + 
  scale_colour_brewer(palette="Set1") +
  my_theme + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))

pdf("Fig3_2hard.pdf")
p
dev.off()

dev.off()
p
