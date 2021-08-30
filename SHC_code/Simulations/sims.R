#' Simulation settings designed for comparing performance of 
#' shc, pvclust, BootClust
#'
#' @param 1 simgroup: set of shc_settings specified below
#' @param 2 sim_iter: the loop index for running multiple instances of same
#'                    simulation to speed up total time - prevents
#'                    multiple instances of same simgroup from overwriting
#'                    each other when running in parallel. Simulations are
#'                    numbered as (sim_iter - 1)*n_reps + i_ter
#' @param 3 n_reps: number of replications of simulation setting to run (e.g. 10)
#' 
#' @author Patrick Kimes

## #############################################################################
## load necessary libraries and functions
## #############################################################################

#setwd("C:/Users/Steve/OneDrive - Imperial College London/MSc/Summer Project/R Code/shc - Copy/sims")

library("methods") #not auto loaded by Rscript
library("R.utils")
#library(parallel) # for detectCores()

#numCores = detectCores()

## parse input parameters
args <- commandArgs(trailingOnly = TRUE)
simgroup <- as.numeric(args[1])
sim_iter <- as.numeric(args[2])
n_reps <- as.numeric(args[3])

#sim_iter = 1
#n_reps = 1

## load testing package
#devtools::load_all("../pkgSHC/")
devtools::load_all("pkgSHC/")

## functions 
source('shc_simulator.R')
source('shc_settings.R')

## nest into linkage and n_bs folder
out_dir <- file.path("outputs")
dir.create(out_dir, showWarnings=FALSE)

## all simulations _FIXXXX
allsims <- c(11001:11045,
             21001:21075, 22001:22075,
             31001:31075, 32001:32075, 33001:33075, 34001:34075, 
             41001:41075,
             51001:51075)


## #############################################################################
## beginning of simulation loops
## #############################################################################
#for (simgroup in 46:195){
## determine which simulations to use
ipart <- allsims[simgroup]

## fix seed value
set.seed(12202015 + sim_iter + 1000 * simgroup)
ipart_dir <- file.path(out_dir, paste0('ipart', ipart))
dir.create(ipart_dir, showWarnings=FALSE)

pZ_shc2CI_1 <- c()
pQ_shc2CI_1 <- c()
time_shc2CI_1 <- c()

pZ_shc2LI_1 <- c()
pQ_shc2LI_1 <- c()
time_shc2LI_1 <- c()

pZ_shc2CI_2 <- c()
pQ_shc2CI_2 <- c()
time_shc2CI_2 <- c()

pZ_shc2LI_2 <- c()
pQ_shc2LI_2 <- c()
time_shc2LI_2 <- c()

cutoff_shc <- c()

true_labs <- c()
hc_list <- list()

for (isim in 1:n_reps) {
    ## simulate data as (d x n) matrix
    output <- shc_settings(ipart)
    data <- output$data
    
    ## save true labels and hclust result
    true_labs <- cbind(true_labs, output$labels)
    
    if (nrow(data) <= 20) {
        ## apply shc() with icovest = 2 (sample)
        start <- proc.time()
        shcOutput1_1 <- shc(t(data), metric="euclidean", linkage="ward.D2",
                            icovest=2, rcpp=FALSE, alpha=0.05,
                            ci=c("2CI"), null_alg=c("hclust"))
        time1 <- proc.time() - start
        
        start <- proc.time()
        shcOutput1_2 <- shc(t(data), metric="euclidean", linkage="ward.D2",
                            icovest=2, rcpp=FALSE, alpha=0.05,
                            ci=c("linkage"), null_alg=c("hclust"))
        time2 <- proc.time() - start
        
        ## save hclust result
        hc_list[[isim]] <- shcOutput1_1$hc_dat ## hc_dat equiv across runs
        
        ## 2-means CI based testing
        pZ_shc2CI_1 <- cbind(pZ_shc2CI_1, shcOutput1_1$p_norm)
        pQ_shc2CI_1 <- cbind(pQ_shc2CI_1, shcOutput1_1$p_emp)
        ## linkage value based testing
        pZ_shc2LI_1 <- cbind(pZ_shc2LI_1, shcOutput1_2$p_norm)
        pQ_shc2LI_1 <- cbind(pQ_shc2LI_1, shcOutput1_2$p_emp)
        ## general testing output
        cutoff_shc <- cbind(cutoff_shc, fwer_cutoff(shcOutput1_1, 0.05))
        time_shc2CI_1 <- rbind(time_shc2CI_1, time1)
        time_shc2LI_1 <- rbind(time_shc2LI_1, time2)
        
    } else {
        start <- proc.time()
        shcOutput2_1 <- shc(t(data), metric="euclidean", linkage="ward.D2",
                            icovest=1, rcpp=FALSE, alpha=0.05,
                            ci=c("2CI"), null_alg=c("hclust"))
        time1 <- proc.time() - start
        
        start <- proc.time()
        shcOutput2_2 <- shc(t(data), metric="euclidean", linkage="ward.D2",
                            icovest=1, rcpp=FALSE, alpha=0.05,
                            ci=c("linkage"), null_alg=c("hclust"))
        time2 <- proc.time() - start
        
        ## save hclust result
        hc_list[[isim]] <- shcOutput2_1$hc_dat ## hc_dat equiv across runs
        
        ## 2-means CI based testing
        pZ_shc2CI_2 <- cbind(pZ_shc2CI_2, shcOutput2_1$p_norm)
        pQ_shc2CI_2 <- cbind(pQ_shc2CI_2, shcOutput2_1$p_emp)
        ## linkage value based testing
        pZ_shc2LI_2 <- cbind(pZ_shc2LI_2, shcOutput2_2$p_norm)
        pQ_shc2LI_2 <- cbind(pQ_shc2LI_2, shcOutput2_2$p_emp)
        ## general testing output
        cutoff_shc <- cbind(cutoff_shc, fwer_cutoff(shcOutput2_1, 0.05))
        time_shc2CI_2 <- rbind(time_shc2CI_2, time1)
        time_shc2LI_2 <- rbind(time_shc2LI_2, time2)
    }
}

## helper function to write results to table
Write2File <- function(ipart, table, sim_iter) {
    fname <- paste0('sim', ipart, '_', deparse(substitute(table)),
                    '_p', sim_iter, '.txt')
    write.table(table,
                file.path(ipart_dir, fname), quote=FALSE,
                sep='\t', row.names=FALSE, col.names=FALSE)
}

## write generating labels to file
Write2File(ipart, true_labs, sim_iter)
save(hc_list, file=file.path(ipart_dir,
                             paste0('sim', ipart, '_hc_list_p', sim_iter, '.rdata')))


Write2File(ipart, cutoff_shc, sim_iter)

## write SHC simulation output to file  
if (nrow(data) <= 20) {
    Write2File(ipart, pZ_shc2CI_1, sim_iter)
    Write2File(ipart, pQ_shc2CI_1, sim_iter)
    Write2File(ipart, pZ_shc2LI_1, sim_iter)
    Write2File(ipart, pQ_shc2LI_1, sim_iter)
    Write2File(ipart, time_shc2CI_1, sim_iter)
    Write2File(ipart, time_shc2LI_1, sim_iter)
} else {
    Write2File(ipart, pZ_shc2CI_2, sim_iter)
    Write2File(ipart, pQ_shc2CI_2, sim_iter)
    Write2File(ipart, pZ_shc2LI_2, sim_iter)
    Write2File(ipart, pQ_shc2LI_2, sim_iter)
    Write2File(ipart, time_shc2CI_2, sim_iter)
    Write2File(ipart, time_shc2LI_2, sim_iter)
}

## let me know that a simulation setting has ended
cat(paste("finished with ipart ", ipart, "\n"))  
#}

