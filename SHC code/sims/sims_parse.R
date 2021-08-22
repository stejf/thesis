#' code to parse output of sims.R for final SHC manuscript tables
#'
#' @param 1 dirname: directory containing output subdirectories of the form
#'          ipart11001/, ipart11002/, ..
#' @param 2 K: simulation group
#' @param 3 n_files: range of sim_iter values used when running job
#' @param 4 n_reps: number of replications of simulation setting per file
#' @param 5 write_sfx: suffix to be used for output table
#' 
#' @author Patrick Kimes


## #############################################################################
## load necessary libraries and functions
## #############################################################################

library("methods") ## not auto loaded by Rscript
library("Hmisc")
library("reshape2")
library("WGCNA") ## for computing Rand Index

## determine currect working directory
#setwd("C:/Users/Steve/OneDrive - Imperial College London/MSc/Summer Project/R Code/shc - Copy/sims")
cwd <- getwd()

args <- commandArgs(trailingOnly = TRUE)
result_dir <- args[1]
K <- as.numeric(args[2])
n_files <- as.numeric(args[3])
n_reps <- as.numeric(args[4])
write_sfx <- args[5]
#result_dir <- "output_mine"
#K <- 2
#n_files <- 1
#n_reps <- 1
#write_sfx <- "_test"

n_treps <- n_files * n_reps

## load helper functions
source("shc_settings.R")
source("shc_simulator.R")
source("shc_helpers.R")



## #############################################################################
## preliminary variables for reading in simulation output
## #############################################################################

simset_list <- list(11001:11045,
                    c(21001:21075, 22001:22075),
                    c(31001:31075, 32001:32075, 33001:33075, 34001:34075),
                    41001:41075,
                    51001:51075)



## pull list of simulation settings
simset <- simset_list[[K]]

## specify names of read-in files (1 = covest by sample)
filenames1 <- c("pZ_shc2LI_1", "time_shc2LI_1",
                "pZ_shc2CI_1", "time_shc2CI_1")

## specify names of read-in files (2 = covest by soft thresholding)
filenames2 <- c("pZ_shc2LI_2", "time_shc2LI_2", 
                "pZ_shc2CI_2", "time_shc2CI_2")

## summary table of all results
summary_table <- data.frame(ipart = numeric(), 
                            n = numeric(),
                            p = numeric(),
                            v = numeric(),
                            delta = numeric(),
                            shape = character(),
                            df = numeric(),
                            ##
                            SHC2zlink = vector(),
                            SHC2zCI = vector(),
                            ##
                            SHC2link_time = numeric(),
                            SHC2CI_time = numeric(), 
                            ## 
                            SHC2zlink_avg = vector(),
                            SHC2zCI_avg = vector(),
                            ## 
                            aRI = vector(),
                            SHC2zlink_aRI = vector(),
                            SHC2zCI_aRI = vector(),
                            stringsAsFactors = FALSE)



## #############################################################################
## read in result of each simulation
## #############################################################################
for (idx in 1:length(simset)) {
    
    ## determine simulation number
    ipart <- simset[idx]
    print(paste("running:", ipart))

    ## add simulation information to table
    siminfo <- shc_settings(ipart, TRUE)
    summary_table[idx, 'ipart'] <- ipart
    summary_table[idx, 'n'] <- siminfo$ni[[1]]*K
    summary_table[idx, 'p'] <- siminfo$d
    summary_table[idx, 'v'] <- siminfo$vsd[1]^2
    summary_table[idx, 'delta'] <- siminfo$delta
    if (length(siminfo$ni) > 1) {
        summary_table[idx, 'shape'] <- paste0(siminfo$shape, "-", paste0(siminfo$ni[[2]], collapse="/"))
    } else {
        summary_table[idx, 'shape'] <- siminfo$shape
    }
    summary_table[idx, 'df'] <- siminfo$df
    
    ## for high-dim use soft threshold, low-dim (if no. of var. p <=20) use sample cov
    filenames <- if (siminfo$d <= 20) { filenames1 } else { filenames2 }

    
    ## gather p-values for SHC methods
    output <- get_results(ipart,
                          filenames=c(filenames[c(1, 3)], "cutoff_shc"),
                          cutoffs=TRUE, nfiles=n_files, pathbase=result_dir) 
    allresults <- output$allresults
    listresults <- output$listresults
    allcutoffs <- output$allcutoffs

    ## skip if simulation output doesn't exist
    if (is.null(allresults)) {
        next
    }
    
    ## determine number of significant calls at root node
    pred_nclusters <- matrix(ncol=ncol(allcutoffs), 
                             nrow=length(listresults))
    for (ifile in 1:length(listresults)) {
        pred_nclusters[ifile, ] <- apply(listresults[[ifile]] <= allcutoffs,
                                         2, function(x) which.min(rev(x)))
    }

    
    ## handle different outputs depending on whether K >< 2
    if (K == 1) {
        call_sig <- rowSums(pred_nclusters > 1)
        summary_table[idx, c('SHC2zlink', 'SHC2zCI')] <- call_sig
    } else {
        ## count number of times predicted number of clusters is correct
        call_k <- rowSums(pred_nclusters == K)
        summary_table[idx, c('SHC2zlink', 'SHC2zCI')] <- call_k
    }

    ## compute the mean predicted number of clusters
    call_kavg <- round(rowMeans(pred_nclusters), 2)
    summary_table[idx, c('SHC2zlink_avg', 'SHC2zCI_avg')] <- call_kavg

    ## gather timing information
    output <- get_times(ipart, filenames[grepl('time_shc2LI', filenames)], nfiles=n_files,
                        pathbase=result_dir)
    summary_table[idx, 'SHC2link_time'] <- round(median(output[, 3]), 2)
    output <- get_times(ipart, filenames[grepl('time_shc2CI', filenames)], nfiles=n_files,
                        pathbase=result_dir)
    summary_table[idx, 'SHC2CI_time'] <- round(median(output[, 3]), 2)
    
    
    ## read in table of true labels
    true_labs <- c()
    for (i in 1:n_files) {
        temp <- read.table(file.path(result_dir, paste0("ipart", ipart),
                                     paste0("sim", ipart, "_true_labs_p", i, ".txt")))
        true_labs <- cbind(true_labs, as.matrix(temp))
    }

    ## read in predicted labels for specific K using dendrogram data object
    pred_labs <- c()
    pred_SHCL <- c()
    pred_SHC2 <- c()
    for (i in 1:n_files) {
        load(file.path(result_dir, paste0("ipart", ipart),
                       paste0("sim", ipart, "_hc_list_p", i, ".rdata")))

        temp <- sapply(hc_list, cutree, K)
        pred_labs <- cbind(pred_labs, temp)

        temp <- mapply(function(x, nc) {cutree(x, nc)}, hc_list,
                       pred_nclusters[1, ((i-1)*n_reps)+(1:n_reps)])
        pred_SHCL <- cbind(pred_SHCL, temp)
        
        temp <- mapply(function(x, nc) {cutree(x, nc)}, hc_list,
                       pred_nclusters[2, ((i-1)*n_reps)+(1:n_reps)])
        pred_SHC2 <- cbind(pred_SHC2, temp)
    }

    ## calculate Adjusted Rand Index (ARI) for all simulations
    ari <- mapply(function(x, y) randIndex(table(x, y)),
                  as.data.frame(true_labs),
                  as.data.frame(pred_labs))
    #to understand better:
    #table(as.data.frame(true_labs)[, 1], as.data.frame(pred_labs)[, 1])

    ## calculate mean ARI for all sims based on cutree on hierarchical solution
    summary_table[idx, 'aRI'] <- round(mean(ari), 2)

    ## calculate mean ARI for SHC_L
    SHCL_ari <- mapply(function(x, y) randIndex(table(x, y)),
                       as.data.frame(true_labs),
                       as.data.frame(pred_SHCL))
    summary_table[idx, 'SHC2zlink_aRI'] <- round(mean(SHCL_ari), 2)


    ## calculate mean ARI for SHC_2
    SHC2_ari <- mapply(function(x, y) randIndex(table(x, y)),
                       as.data.frame(true_labs),
                       as.data.frame(pred_SHC2))
    summary_table[idx, 'SHC2zCI_aRI'] <- round(mean(SHC2_ari), 2)
}

## table to csv
write.csv(summary_table, file=file.path(cwd, paste0("K", K, "_summary_table", write_sfx, ".csv")))


## clean up table
if (K == 1) {
    print_table <- summary_table[, c("ipart", "n", "p", "v", "delta", "shape", "df",
                                     "SHC2zlink", "SHC2zCI", 
                                     "SHC2link_time", "SHC2CI_time")]
    ## table to csv
    write.csv(print_table, file=file.path(cwd, paste0("cleaned-K", K, "_summary_table", write_sfx, ".csv")))
    
} else {
    summary_table$pr_SHC2zlink <- mapply(function(x, y) { paste0(x, " (", y, ")") },
                                         summary_table$SHC2zlink, summary_table$SHC2zlink_avg)
    summary_table$pr_SHC2zCI <- mapply(function(x, y) { paste0(x, " (", y, ")") },
                                       summary_table$SHC2zCI, summary_table$SHC2zCI_avg)

    print_table <- summary_table[, c("ipart", "n", "p", "v", "delta", "shape", "df",
                                     "pr_SHC2zlink", "pr_SHC2zCI",
                                     "SHC2link_time", "SHC2CI_time", 
                                     "aRI", "SHC2zlink_aRI", "SHC2zCI_aRI")]
    ## table to csv
    write.csv(print_table, file=file.path(cwd, paste0("cleaned-K", K, "_summary_table", write_sfx, ".csv")))
}




