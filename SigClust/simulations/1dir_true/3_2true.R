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
n = 10
d = 12
mu = rep(0, d)
a = c(0, 10, 20, 25)
D = diag(d)
v = 35 
D[1, 1] = v
D[1, ]
nsim = 1000

CI_original1 = foreach(b=1:100, .packages='MASS', .combine='c') %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D)
  dat[(n/2+1):(n), 1] = dat[(n/2+1):(n), 1] + a[1]
  km = kmeans(dat, 2, nstart=20)
  CI_original = km$tot.withinss / km$totss
  return(CI_original)
}

D_star = D
D_star[1, 1] = D[1, 1] + 0.25*a[1]^2

simulations1 = foreach (j=101:1100, .combine='c', .packages='MASS') %dopar%{
  set.seed(1+j)
  dat_sim = mvrnorm(n=n, mu=mu, Sigma=D_star)
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

CI_original2 = foreach(b=1:100, .packages='MASS', .combine='c') %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D)
  dat[(n/2+1):(n), 1] = dat[(n/2+1):(n), 1] + a[2]
  km = kmeans(dat, 2, nstart=20)
  CI_original = km$tot.withinss / km$totss
  return(CI_original)
}

D_star = D
D_star[1, 1] = D[1, 1] + 0.25*a[2]^2

simulations2 = foreach (j=101:1100, .combine='c', .packages='MASS') %dopar%{
  set.seed(1+j)
  dat_sim = mvrnorm(n=n, mu=mu, Sigma=D_star)
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

CI_original3 = foreach(b=1:100, .packages='MASS', .combine='c') %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D)
  dat[(n/2+1):(n), 1] = dat[(n/2+1):(n), 1] + a[3]
  km = kmeans(dat, 2, nstart=20)
  CI_original = km$tot.withinss / km$totss
  return(CI_original)
}

D_star = D
D_star[1, 1] = D[1, 1] + 0.25*a[3]^2

simulations3 = foreach (j=101:1100, .combine='c', .packages='MASS') %dopar%{
  set.seed(1+j)
  dat_sim = mvrnorm(n=n, mu=mu, Sigma=D_star)
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

CI_original4 = foreach(b=1:100, .packages='MASS', .combine='c') %dopar%{
  set.seed(1+b)
  dat = mvrnorm(n=n, mu=mu, Sigma=D)
  dat[(n/2+1):(n), 1] = dat[(n/2+1):(n), 1] + a[4]
  km = kmeans(dat, 2, nstart=20)
  CI_original = km$tot.withinss / km$totss
  return(CI_original)
}

D_star = D
D_star[1, 1] = D[1, 1] + 0.25*a[4]^2

simulations4 = foreach (j=101:1100, .combine='c', .packages='MASS') %dopar%{
  set.seed(1+j)
  dat_sim = mvrnorm(n=n, mu=mu, Sigma=D_star)
  km_sim = kmeans(dat_sim, 2, nstart=20)
  CI_sim = km_sim$tot.withinss / km_sim$totss
  return(CI_sim)
}

CI_sim4 = simulations4
pvalue4 = c()
for (i in 1:100){
  pvalue4[i] = (1/nsim)*sum(as.numeric(CI_sim4 < CI_original4[i]))
}

save(pvalue1, pvalue2, pvalue3, pvalue4, file="pvalues.RData")
###-------------------------------------------------------------------------- 
#Plots
rm(list=ls())
setwd('C:/Users/Steve/OneDrive - Imperial College London/MSc/Summer Project/R Code/Huang d_300/1dir_true')
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

# mean p-values
mean(pvalue1)
mean(pvalue2)
mean(pvalue3)
mean(pvalue4)

# number of p-values less than 0.05
sum(pvalue1 < 0.05)
sum(pvalue2 < 0.05)
sum(pvalue3 < 0.05)
sum(pvalue4 < 0.05)

# number of p-values less than 0.1
sum(pvalue1 < 0.1)
sum(pvalue2 < 0.1)
sum(pvalue3 < 0.1)
sum(pvalue4 < 0.1)

pvals = as.data.frame(cbind(pvalue1, pvalue2, pvalue3, pvalue4))
names(pvals) = a
pvals = melt(pvals, measure.vars=as.character(a), variable.name="a", 
             value.name="pvalue")

p = ggplot() + stat_ecdf(dat=pvals, aes(x=pvalue, linetype=a, colour=a)) + 
  xlim(0, 1) + geom_vline(xintercept=0.05, color="gray45") +
  labs(x="p-value", y="Empirical Distribution") + 
  ggtitle("True Method (Mixture 2 Gauss. - 1 Direction)") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate(geom = "text", label = bquote(alpha == .(0.05)), x = 0.05, y = 1.1, 
           size=5) + 
  scale_colour_brewer(palette="Set1") +
  my_theme + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))

pdf("Fig3_2true.pdf")
p
dev.off()

dev.off()
p
