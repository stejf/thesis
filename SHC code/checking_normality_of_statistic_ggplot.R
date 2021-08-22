rm(list=ls())
setwd("C:/Users/Steve/OneDrive - Imperial College London/MSc/Summer Project/R Code/shc - Copy/sims")
devtools::load_all("../pkgSHC/")
source('shc_simulator.R')
source('shc_settings.R')

my_theme = theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12))

set.seed(2)
# pick any dataset by changing ipart - refer to shc_settings 
# for all iparts
ipart <- 31005 # just picking a dataset - different ones were considered
output <- shc_settings(ipart)
data <- output$data
x <- t(data)
metric = "euclidean"
linkage = "ward.D2"
l = 2
alpha = 0.05
icovest = 1
bkgd_pca = TRUE
n_min = 10
rcpp = FALSE
ci = "2CI"
null_alg = "hclust"
ci_idx = 1
ci_emp = FALSE
n <- nrow(x)
p <- ncol(x)
n_ci <- length(ci)
x_clust <- .initcluster(x, n, p, metric, linkage, l, 
                        n_ci, ci, rcpp)
ci_dat <- x_clust$ci_dat
hc_dat <- x_clust$hc_dat
idx_hc <- x_clust$idx_hc

##p-values for all <=(n-1) tests
p_emp <- matrix(2, nrow=n-1, ncol=n_ci)
p_norm <- matrix(2, nrow=n-1, ncol=n_ci)
colnames(p_emp) <- paste(null_alg, ci, sep="_")
colnames(p_norm) <- paste(null_alg, ci, sep="_")

##null covariance parameters for all <=(n-1) tests
eigval_dat <- matrix(-1, nrow=n-1, ncol=p)
eigval_sim <- matrix(-1, nrow=n-1, ncol=p)
backvar <- rep(-1, n-1)

##############################################################################
# First we compare the distributions when n_sim=100 vs when n_sim=1000. 
# The plot is saved in variable 'result' with the former density displayed in
# black, and the latter displayed in red.

# CHOOSE NODE k - 100 simulations ---------------------------------------------
k= n-1 # different nodes were considered here
n_sim = 100
ci_sim <- array(-1, dim=c(n-1, n_sim, n_ci))

##indices for subtree
idx_sub <- unlist(idx_hc[k, ])
n_sub <- length(idx_sub)

##estimate null Gaussian
xk_null <- null_eigval(x[idx_sub, ], n_sub, p, icovest, bkgd_pca)

##prevent messages from looped application of clustering
## messages will be thrown once at initial clustering
suppressMessages(
  ##simulate null datasets
  for (i in 1:n_sim) {
    xsim <- .simnull(xk_null$eigval_sim, n_sub, p)
    ci_sim[k, i, ] <- .calcCI_shc(xsim, p, metric, linkage, l, 
                                  n_ci, ci, null_alg, rcpp)
  }
)


A = ci_sim[k, , ]


## fit density to cluster indices
den_pv_A <- density(A)
den_df_A <- data.frame(x=den_pv_A$x, y=den_pv_A$y)

## determine range of density plot
xmin_A <- min(den_df_A$x)
ymax_A <- max(den_df_A$y)

mindex_A <- mean(A)
sindex_A <- sd(A)

gpA <- ggplot() +
  geom_path(aes(x=x, y=y), data=den_df_A, color="red", size=1) +
  geom_vline(xintercept=ci_dat[k], color="#4daf4a", size=1) + 
  ylab("density") + xlab("Cluster Index") + 
  ggtitle(paste0("SHC Results for node ", n-k)) +
  scale_x_continuous(expand=c(0.01, 0)) +
  scale_y_continuous(expand=c(0.01, 0)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  my_theme
print(gpA)
# k- 1000 simulations -----------------------------------------------
n_sim=1000
ci_sim <- array(-1, dim=c(n-1, n_sim, n_ci))

##prevent messages from looped application of clustering
## messages will be thrown once at initial clustering
suppressMessages(
  ##simulate null datasets
  for (i in 1:n_sim) {
    xsim <- .simnull(xk_null$eigval_sim, n_sub, p)
    ci_sim[k, i, ] <- .calcCI_shc(xsim, p, metric, linkage, l, 
                                  n_ci, ci, null_alg, rcpp)
  }
)
B = ci_sim[k, , ]

## fit density to cluster indices
den_pv_B <- density(B)
den_df_B <- data.frame(x=den_pv_B$x, y=den_pv_B$y)

## determine range of density plot
xmin_B <- min(den_df_B$x)
ymax_B <- max(den_df_B$y)

mindex_B <- mean(B)
sindex_B <- sd(B)

result = gpA + geom_path(aes(x=x, y=y), data=den_df_B, color="black", size=1)
print(result)
pdf("densities2.pdf")
result
dev.off()
# recall that since we are using 2CI, then the green line would need to be to
# the left of the densities for us to have a chance at rejecting
# - so if we have a clear case of just one cluster the green line will
# be to the right

##############################################################################
# Next we compare both distributions of n_sim=100 and n_sim=1000, with the 
# approximate normal distribution that would be used instead.
# These two plots are saved in 'gpAA' and 'gpBB'

######## n_sim=100
n_sim = 100
kci_df_A <- data.frame(x=A, y=runif(length(A), ymax_A*1/4, ymax_A*3/4))

##compute p-values
m_idx_A <- colMeans(as.matrix(A))
s_idx_A <- apply(as.matrix(A), 2, sd)
p_norm[k, ] <- pnorm(ci_dat[k, ], m_idx_A, s_idx_A)
p_emp[k, ] <- colMeans(as.matrix(A) <= 
                         matrix(ci_dat[k, ], nrow=n_sim, 
                                ncol=n_ci, byrow=TRUE))

gpAA <- ggplot(kci_df_A) +
  geom_point(aes(x=x, y=y), alpha=1/2, size=1, color="#377eb8") +
  geom_path(aes(x=x, y=y), data=den_df_A, color="red", size=1) + 
  stat_function(fun=dnorm, args=list(mean=mindex_A, sd=sindex_A),
                color="red", linetype="dashed", size=1, n=500) +
  geom_vline(xintercept=ci_dat[k], color="#4daf4a", size=1) + 
  ylab("density") + xlab("Cluster Index") + 
  ggtitle(paste0("SHC Results for node ", n-k)) +
  scale_x_continuous(expand=c(0.01, 0)) +
  scale_y_continuous(expand=c(0.01, 0)) + 
  annotate(geom="text", x=min(xmin_A, ci_dat[k]), y=ymax_A*.98, vjust=1, hjust=0,
               label=paste0(" p-value (Q) = ", round(p_emp[k, ], 3), "\n",
                            " p-value (Z) = ", round(p_norm[k, ], 3))) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  my_theme
print(gpAA)
pdf("densities3.pdf")
gpAA
dev.off()

######## n_sim=1000
n_sim = 1000
kci_df_B <- data.frame(x=B, y=runif(length(B), ymax_B*1/4, ymax_B*3/4))
##compute p-values
m_idx_B <- colMeans(as.matrix(B))
s_idx_B <- apply(as.matrix(B), 2, sd)
p_norm[k, ] <- pnorm(ci_dat[k, ], m_idx_B, s_idx_B)
p_emp[k, ] <- colMeans(as.matrix(B) <= 
                         matrix(ci_dat[k, ], nrow=n_sim, 
                                ncol=n_ci, byrow=TRUE))

gpBB <- ggplot(kci_df_B) +
  geom_point(aes(x=x, y=y), alpha=1/2, size=1, color="#377eb8") +
  geom_path(aes(x=x, y=y), data=den_df_B, color="black", size=1) + 
  stat_function(fun=dnorm, args=list(mean=mindex_B, sd=sindex_B),
                color="black", linetype="dashed", size=1, n=500) +
  geom_vline(xintercept=ci_dat[k], color="#4daf4a", size=1) + 
  ylab("density") + xlab("Cluster Index") + 
  ggtitle(paste0("SHC Results for node ", n-k)) +
  scale_x_continuous(expand=c(0.01, 0)) +
  scale_y_continuous(expand=c(0.01, 0)) + annotate(geom="text", x=min(xmin_B, ci_dat[k]), y=ymax_B*.98, vjust=1, hjust=0,
                label=paste0(" p-value (Q) = ", round(p_emp[k, ], 3), "\n",
                             " p-value (Z) = ", round(p_norm[k, ], 3))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  my_theme
print(gpBB)
pdf("densities1.pdf")
gpBB
dev.off()


##############################################################################
# we test for normality when n_sim=100 (case A) and when n_sim=1000 (case B)
a=ks.test(A, "pnorm", mean(A), sd(A))
b=ks.test(B, "pnorm", mean(B), sd(B))
a
b

# we obtain qqplots for n_sim=100 (case A) and when n_sim=1000 (case B)
qqnorm(A)
qqnorm(B)
