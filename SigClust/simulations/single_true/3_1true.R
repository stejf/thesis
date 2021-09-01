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
nsim = 1000

CI_original1 = foreach(b=1:100, .packages='MASS', .combine='c') %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D)
  km = kmeans(dat, 2, nstart=20)
  CI_original = km$tot.withinss / km$totss
  return(CI_original)
}

simulations1 = foreach (j=101:1100, .combine='c', .packages='MASS') %dopar%{
  set.seed(1+j)
  dat_sim = mvrnorm(n=n, mu=mu, Sigma=D)
  km_sim = kmeans(dat_sim, 2, nstart=20)
  CI_sim = km_sim$tot.withinss / km_sim$totss
  return(CI_sim)
}

CI_sim1 = simulations1
pvalue1 = c()
for (i in 1:100){
  pvalue1[i] = (1/nsim)*sum(as.numeric(CI_sim1 < CI_original1[i]))
}

###-------------------------------------------------------------------------- 
D[1, 1] = v[2]
D[1, ]

CI_original2 = foreach(b=1:100, .packages='MASS', .combine='c') %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D)
  km = kmeans(dat, 2, nstart=20)
  CI_original = km$tot.withinss / km$totss
  return(CI_original)
}

simulations2 = foreach (j=101:1100, .combine='c', .packages='MASS') %dopar%{
  set.seed(1+j)
  dat_sim = mvrnorm(n=n, mu=mu, Sigma=D)
  km_sim = kmeans(dat_sim, 2, nstart=20)
  CI_sim = km_sim$tot.withinss / km_sim$totss
  return(CI_sim)
}

CI_sim2 = simulations2
pvalue2 = c()
for (i in 1:100){
  pvalue2[i] = (1/nsim)*sum(as.numeric(CI_sim2 < CI_original2[i]))
}

###-------------------------------------------------------------------------- 
D[1, 1] = v[3]
D[1, ]

CI_original3 = foreach(b=1:100, .packages='MASS', .combine='c') %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D)
  km = kmeans(dat, 2, nstart=20)
  CI_original = km$tot.withinss / km$totss
  return(CI_original)
}

simulations3 = foreach (j=101:1100, .combine='c', .packages='MASS') %dopar%{
  set.seed(1+j)
  dat_sim = mvrnorm(n=n, mu=mu, Sigma=D)
  km_sim = kmeans(dat_sim, 2, nstart=20)
  CI_sim = km_sim$tot.withinss / km_sim$totss
  return(CI_sim)
}

CI_sim3 = simulations3
pvalue3 = c()
for (i in 1:100){
  pvalue3[i] = (1/nsim)*sum(as.numeric(CI_sim3 < CI_original3[i]))
}

###-------------------------------------------------------------------------- 
D[1, 1] = v[4]
D[1, ]

CI_original4 = foreach(b=1:100, .packages='MASS', .combine='c') %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D)
  km = kmeans(dat, 2, nstart=20)
  CI_original = km$tot.withinss / km$totss
  return(CI_original)
}

simulations4 = foreach (j=101:1100, .combine='c', .packages='MASS') %dopar%{
  set.seed(1+j)
  dat_sim = mvrnorm(n=n, mu=mu, Sigma=D)
  km_sim = kmeans(dat_sim, 2, nstart=20)
  CI_sim = km_sim$tot.withinss / km_sim$totss
  return(CI_sim)
}

CI_sim4 = simulations4
pvalue4 = c()
for (i in 1:100){
  pvalue4[i] = (1/nsim)*sum(as.numeric(CI_sim4 < CI_original4[i]))
}

###-------------------------------------------------------------------------- 
D[1, 1] = v[5]
D[1, ]

CI_original5 = foreach(b=1:100, .packages='MASS', .combine='c') %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D)
  km = kmeans(dat, 2, nstart=20)
  CI_original = km$tot.withinss / km$totss
  return(CI_original)
}

simulations5 = foreach (j=101:1100, .combine='c', .packages='MASS') %dopar%{
  set.seed(1+j)
  dat_sim = mvrnorm(n=n, mu=mu, Sigma=D)
  km_sim = kmeans(dat_sim, 2, nstart=20)
  CI_sim = km_sim$tot.withinss / km_sim$totss
  return(CI_sim)
}

CI_sim5 = simulations5
pvalue5 = c()
for (i in 1:100){
  pvalue5[i] = (1/nsim)*sum(as.numeric(CI_sim5 < CI_original5[i]))
}

save(pvalue1, pvalue2, pvalue3, pvalue4, pvalue5, file="pvalues.RData")
###-------------------------------------------------------------------------- 
#Plots
rm(list=ls())
setwd('C:/Users/Steve/OneDrive - Imperial College London/MSc/Summer Project/R Code/Huang d_300/single_true')
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

# mean p-values
mean(pvalue1)
mean(pvalue2)
mean(pvalue3)
mean(pvalue4)
mean(pvalue5)

# number of p-values less than 0.05
sum(pvalue1 < 0.05)
sum(pvalue2 < 0.05)
sum(pvalue3 < 0.05)
sum(pvalue4 < 0.05)
sum(pvalue5 < 0.05)

# number of p-values less than 0.1
sum(pvalue1 < 0.1)
sum(pvalue2 < 0.1)
sum(pvalue3 < 0.1)
sum(pvalue4 < 0.1)
sum(pvalue5 < 0.1)

pvals = as.data.frame(cbind(pvalue1, pvalue2, pvalue3, pvalue4, pvalue5))
names(pvals) = v
pvals = melt(pvals, measure.vars=as.character(v), variable.name="v", 
             value.name="pvalue")

p = ggplot() + stat_ecdf(dat=pvals, aes(x=pvalue, linetype=v, col=v)) + 
  xlim(0, 1) + geom_vline(xintercept=0.05, color="gray45") +
  labs(x="p-value", y="Empirical Distribution") + 
  ggtitle("True Method (Single Gaussian)") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate(geom = "text", label = bquote(alpha == .(0.05)), x = 0.05, y = 1.1, 
           size=5) + 
  scale_colour_brewer(palette="Set1", direction=-1) +
  my_theme + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))

pdf("Fig3_1true.pdf")
p
dev.off()

dev.off()
p


#pvalue = c()
#for(b in 1:4){
#  set.seed(1+b)
#  dat = mvrnorm(n=n, mu=mu, Sigma=D)
#  km = kmeans(dat, 2, nstart=20)
#  CI_original = km$tot.withinss / km$totss
#  simulations = foreach (j=1:1000, .combine='c') %dopar%{  
#    library(MASS) # for mvrnorm() 
#    d = 14
#    mu = rep(0, d)
#    set.seed(1+j)
#    dat_sim = mvrnorm(n=n, mu=mu, Sigma=D)
#    km_sim = kmeans(dat_sim, 2, nstart=20)
#    CI_sim = km_sim$tot.withinss / km_sim$totss
#    return(CI_sim)
#  }
#  nsim = 1000
#  CI_sim = simulations
#  pvalue[b] = (1/nsim)*sum(as.numeric(CI_sim < CI_original))
#  return(pvalue)
#}
#pvalue

#pvalue1 = foreach(b=1:100, .packages=c('foreach', 'doParallel', 'MASS')) %dopar%{
#  set.seed(1+b)
#  dat = mvrnorm(n=n, mu=mu, Sigma=D)
#  km = kmeans(dat, 2, nstart=20)
#  CI_original = km$tot.withinss / km$totss
#  simulations = foreach (j=1:1000, .combine='c', .packages='MASS') %dopar%{  
#    d = 7
#    mu = rep(0, d)
#    set.seed(1+j)
#    dat_sim = mvrnorm(n=n, mu=mu, Sigma=D)
#    km_sim = kmeans(dat_sim, 2, nstart=20)
#    CI_sim = km_sim$tot.withinss / km_sim$totss
#    return(CI_sim)
#  }
#  nsim = 1000
#  CI_sim = simulations
#  pvalue = (1/nsim)*sum(as.numeric(CI_sim < CI_original))
#  return(pvalue)
#}