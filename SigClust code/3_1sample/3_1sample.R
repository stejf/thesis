rm(list=ls())
library(sigclust) # for sigclust()
library(MASS) # for mvrnorm() 
library(tictoc) # for tic() and toc()
library(foreach) # for foreach()
library(parallel) # for detectCores()
library(doParallel) # for registerDoParallel()

numCores = detectCores()
registerDoParallel(numCores)

### CASE OF ONE CLUSTER
###-------------------------------------------------------------------------- 
n = 40
d = 300
mu = rep(0, d)
D = diag(d)
v = c(1, 15, 35, 200, 1000) # different values for first eigenvalue
D[1, 1] = v[1] # v = 1
D[1, ]

pvalue1 = foreach (b=1:100, .packages=c('MASS', 'sigclust')) %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=40 d-dimensional observations
  nsim = 1000
  nrep = 1
  icovest = 2
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}

###-------------------------------------------------------------------------- 
D[1, 1] = v[2]
D[1, ]

pvalue2 = foreach (b=1:100, .packages=c('MASS', 'sigclust')) %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=40 d-dimensional observations
  nsim = 1000
  nrep = 1
  icovest = 2
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}
###-------------------------------------------------------------------------- 
D[1, 1] = v[3]
D[1, ]

pvalue3 = foreach (b=1:100, .packages=c('MASS', 'sigclust')) %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=40 d-dimensional observations
  nsim = 1000
  nrep = 1
  icovest = 2
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}

###-------------------------------------------------------------------------- 
D[1, 1] = v[4]
D[1, ]

pvalue4 = foreach (b=1:100, .packages=c('MASS', 'sigclust')) %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=40 d-dimensional observations
  nsim = 1000
  nrep = 1
  icovest = 2
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}

###-------------------------------------------------------------------------- 
D[1, 1] = v[5]
D[1, ]

pvalue5 = foreach (b=1:100, .packages=c('MASS', 'sigclust')) %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D) # n=40 d-dimensional observations
  nsim = 1000
  nrep = 1
  icovest = 2
  pvalue = sigclust(dat, nsim=nsim, nrep=nrep, labflag=0, icovest=icovest)
  return(pvalue)
}

save(pvalue1, pvalue2, pvalue3, pvalue4, pvalue5, file="pvalues.RData")
###-------------------------------------------------------------------------- 
#Plots
rm(list=ls())
setwd('C:/Users/Steve/OneDrive - Imperial College London/MSc/Summer Project/R Code/Huang d_300/3_1sample')
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

v = c(1, 15, 35, 200, 1000)

sum(pvalue1[[1]]@veigval) # sum of sample eigenvalues
sum(pvalue1[[1]]@vsimeigval) # sum of hard-thresholded eigenvalues
# so hard thresholding method gives a larger denominator of equation 16 as 
# pointed out in Huang et al.

pvalue1[[1]]@veigval[1] == pvalue1[[1]]@vsimeigval[1]
# meanwhile the first eigenvalue is the same for both methods as pointed out


# plots for 1 realisation (1 of the 100 datasets) when v = 1
pdf("v1_all.pdf")
plot(pvalue1[[1]], arg="all")
dev.off()

pdf("v15_all.pdf")
plot(pvalue2[[1]], arg="all")
dev.off()

pdf("v35_all.pdf")
plot(pvalue3[[1]], arg="all")
dev.off()

pdf("v200_all.pdf")
plot(pvalue4[[1]], arg="all")
dev.off()

pdf("v1000_all.pdf")
plot(pvalue5[[1]], arg="all")
dev.off()

pval1 = sapply(1:100, function(i) pvalue1[[i]]@pval)
pval2 = sapply(1:100, function(i) pvalue2[[i]]@pval)
pval3 = sapply(1:100, function(i) pvalue3[[i]]@pval)
pval4 = sapply(1:100, function(i) pvalue4[[i]]@pval)
pval5 = sapply(1:100, function(i) pvalue5[[i]]@pval)

# mean p-values
mean(pval1)
mean(pval2)
mean(pval3)
mean(pval4)
mean(pval5)

# number of p-values less than 0.05
sum(pval1 < 0.05)
sum(pval2 < 0.05)
sum(pval3 < 0.05)
sum(pval4 < 0.05)
sum(pval5 < 0.05)

# number of p-values less than 0.1
sum(pval1 < 0.1)
sum(pval2 < 0.1)
sum(pval3 < 0.1)
sum(pval4 < 0.1)
sum(pval5 < 0.1)

pvals = as.data.frame(cbind(pval1, pval2, pval3, pval4, pval5))
names(pvals) = v
pvals = melt(pvals, measure.vars=as.character(v), variable.name="v", 
             value.name="pvalue")

p = ggplot() + stat_ecdf(dat=pvals, aes(x=pvalue, linetype=v, colour=v)) + 
  xlim(0, 1) + geom_vline(xintercept=0.05, color="gray45") +
  labs(x="p-value", y="Empirical Distribution") + 
  ggtitle("Single Gaussian (sample)") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate(geom = "text", label = bquote(alpha == .(0.05)), x = 0.05, y = 1.1, 
           size=5) + 
  scale_colour_brewer(palette="Set1", direction=-1) +
  my_theme + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))

pdf("Fig3_1sample.pdf")
p
dev.off()

dev.off()
p
