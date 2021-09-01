## #############################################################################
## helper functions for parsing output of shc simulation results
## #############################################################################


#' Parse p-values Results for Simulations
#'
#' @param ipart integer value specifying the simulation index to be parsed
#' @param filenames vector of strings specifying the data output files
#'        to be parsed, if 'cutoffs' is TRUE, then the final entry of the vector
#'        should be the string name of the files containing the FWER cutoff values
#' @param cutoffs logical value whether final entry of 'filenames' specifies the
#'        location of files containing the FWER cutoff values 
#' @param nfiles integer value specifying number of independent simulation loops
#'        run (default = 10)
#' @param pathbase string specifying base directory of simulation output for parsing
#'        (default = '.')
#'
#' @author Patrick Kimes
get_results <- function(ipart, filenames, cutoffs, nfiles = 2, pathbase = '.') {

    flist <- grep(paste0("sim", ipart, "_", filenames[2]),
                  list.files(file.path(pathbase, paste0("ipart", ipart))))

    ## only proceed with parsing if directory contains as many parts as expected
    if (length(flist) != nfiles) {
        return(list(allresults=c(), listresults=c(), allcutoffs=c()))
    }

    allresults <- c()
    listresults <- list()
    for (ifile in 1:(length(filenames) - cutoffs)) {
        ffile <- filenames[ifile]
        full_ffile <- file.path(pathbase, paste0('ipart', ipart),
                                paste0('sim', ipart, '_', ffile))
        results <- c()
        for (igroup in 1:nfiles) {
            iresults <- read.table(paste0(full_ffile, '_p', igroup, '.txt'), 
                                   header=FALSE) # (nsamp-1) x nsims
            iresults <- as.matrix(iresults)
            results <- cbind(results, iresults)
        }
        colnames(results) <- paste0('sim', 1:ncol(results))
        
        listresults[[ifile]] <- results
        results <- data.frame(results, 'node'=1:nrow(results))
        results <- melt(results, id.vars='node')
        results <- cbind(results, 'method'=ffile)
        allresults <- rbind(allresults, results)
    }
    
    ## gather FWER control cutoff information
    allcutoffs <- c()
    if (cutoffs) {
        ffile <- tail(filenames, n=1)
        full_ffile <- file.path(pathbase, paste0('ipart', ipart),
                                paste0('sim', ipart, '_', ffile))
        cutoffs <- c()
        for (igroup in 1:nfiles) {
            cutoffs <- read.table(paste0(full_ffile, '_p', igroup, '.txt'),
                                  header=FALSE)
            cutoffs <- as.matrix(cutoffs)
            allcutoffs <- cbind(allcutoffs, cutoffs)
        }
    }
    
    ## return results
    list(allresults=allresults,
         listresults=listresults,
         allcutoffs=allcutoffs)
}



#' Parse Timing Results for Simulations
#'
#' @param ipart integer value specifying the simulation index to be parsed
#' @param filename string specifying the data output files to be parsed
#' @param nfiles integer value specifying number of independent simulation loops
#'        run (default = 10)
#' @param pathbase string specifying base directory of simulation output for parsing
#'        (default = '.')
#'
#' @author Patrick Kimes
get_times <- function(ipart, filename, nfiles, pathbase = '.') {  
    full_ffile <- file.path(pathbase, paste0('ipart', ipart), 
                            paste0('sim', ipart,  '_', filename))
    times <- c()
    for (igroup in 1:nfiles) {
        itimes <- read.table(paste0(full_ffile, '_p', igroup, '.txt'), 
                             header=FALSE) # (nsamp-1) x nsims
        itimes <- as.matrix(itimes)
        times <- rbind(times, itimes)
    }
    times
}


