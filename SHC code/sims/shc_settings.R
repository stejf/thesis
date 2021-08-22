#' @title Simulation Settings for SHC Testing
#' 
#' note: vspikes are on variance scale
#'
#' @param ipart integer specifying simulation setting of interest
#' @param params logical specifying whether to return parameter values
#'        instead of an actual dataset
#'
#' @details
#' Settings (dim=1000 for all cases):
#' 
#' K=1 simulation settings
#'   11001:11045 : spikes=1  [5v x 3df x 3ni]
#' 
#'
#' K=2 simulation settings [line - mean signal in first direction]
#'   21001:21075 : spikes=1  [5v x 5d x 3df x 1ni]
#'
#' K=2 simulation settings [line - mean signal in second direction]
#'   22001:22075 : spikes=1  [5v x 5d x 3df x 1ni]
#'
#' 
#' K=3 simulation settings [line - mean signal in first direction]
#'   31001:31075 : spikes=1  [5v x 5d x 3df x 1ni]
#' 
#' K=3 simulation settings [line - mean signal in second direction]
#'   32001:32075 : spikes=1  [5v x 5d x 3df x 1ni]
#' 
#' K=3 simulation settings [triangle + jitter - mean signal in first 2 directions]
#'   33001:33075 : spikes=1  [5v x 5d x 3df x 1ni]
#' 
#' K=3 simulation settings [triangle + jitter - mean signal in 2nd and 3rd directions]
#'   34001:34075 : spikes=1  [5v x 5d x 3df x 1ni]
#' 
#' 
#' K=4 simulation settings [tetrahedron + jitter]
#'   41001:41075 : spikes=1  [5v x 5d x 3df x 1ni]
#'
#'
#' K=5 simulation settings [random in K-1 dim sphere]
#'   51001:51075 : spikes=1  [5v x 5d x 3df x 1ni] - 2,6,10,14,18
#' 
#' @author Patrick Kimes
shc_settings <- function(ipart, params = FALSE) {

    ## #########################################################################
    ## simulation settings for K=1
    ## #########################################################################
    if (ipart %in% 11000:11999) {
        if (ipart %in% 11001:11045) {
            ## #################################################################
            ## dim = 1000, spikes = 1 
            ## #################################################################
            shape <- "-"
            base <- 11000
            q <- 1
            d <- 700
            delta <- 0
        } else {
            stop("invalid ipart")
        }

        ## shared parameters for all K = 1 simulation settings
        k <- 1
        vspikes <- c(1, 10, 50, 100, 1000)
        dfs <- c(-1, 6, 3)
        nis <- c(50, 100, 200)
        
        ispike <- ((ipart - base - 1) %% 5) + 1
        idf <- (floor((ipart - base - 1) / 5) %% 3) + 1
        ini <- ceiling((ipart - base) / (5*3))

        spike <- vspikes[ispike]
        df <- dfs[idf]
        ni <- nis[ini]

        m1 <- rep(0, d)
        mui <- rbind(m1)
        vsd <- c(rep(sqrt(spike), q), rep(1, d - q))
    }
    
    
    ## #########################################################################
    ## simulation settings for K=2 [line - mean signal in first direction] 
    ## #########################################################################
    if (ipart %in% 21000:21999) {
        if (ipart %in% 21001:21075) {
            ## #################################################################
            ## dim = 1000, spikes = 1
            ## #################################################################
            base <- 21000
            q <- 1
            d <- 700
        } else {
            stop("invalid ipart")
        }

        shape <- "line - mean signal in 1st direction"
        k <- 2
        vspikes <- c(1, 10, 50, 100, 1000)
        dfs <- c(-1, 6, 3)
        ni <- 50
        deltas <- seq(2, 10, by=2)
        
        ispike <- (floor((ipart - base - 1) / 15) %% 5) + 1
        idelta <- ((ipart - base - 1) %% 5) + 1
        idf <- (floor((ipart - base - 1) / 5) %% 3) + 1

        spike <- vspikes[ispike]
        delta <- deltas[idelta]
        df <- dfs[idf]

        m1 <- c(-delta/2, rep(0, d-1))
        m2 <- c( delta/2, rep(0, d-1))
        mui <- rbind(m1, m2)
        vsd <- c(rep(sqrt(spike), q), rep(1, d - q))
    }
    
    ## #########################################################################
    ## simulation settings for K=2 [line - mean signal in second direction]
    ## #########################################################################
    if (ipart %in% 22000:22999) {
        if (ipart %in% 22001:22075) {
            ## #################################################################
            ## dim = 1000, spikes = 1
            ## #################################################################
            base <- 22000
            q <- 1
            d <- 700
        } else {
            stop("invalid ipart")
        }
        
        shape <- "line - mean signal in 2nd direction"
        k <- 2
        vspikes <- c(1, 10, 50, 100, 1000)
        dfs <- c(-1, 6, 3)
        ni <- 50
        deltas <- seq(2, 10, by=2)
        
        ispike <- (floor((ipart - base - 1) / 15) %% 5) + 1
        idelta <- ((ipart - base - 1) %% 5) + 1
        idf <- (floor((ipart - base - 1) / 5) %% 3) + 1

        spike <- vspikes[ispike]
        delta <- deltas[idelta]
        df <- dfs[idf]

        m1 <- c(0, -delta/2, rep(0, d-2))
        m2 <- c(0, delta/2, rep(0, d-2))
        mui <- rbind(m1, m2)
        vsd <- c(rep(sqrt(spike), q), rep(1, d - q))
    }


    ## #########################################################################
    ## K=3 simulation settings [line - mean signal in first direction]
    ## #########################################################################
    if (ipart %in% 31000:31999) {
        if (ipart %in% 31001:31075) {
            ## #################################################################
            ## dim = 1000, spikes = 1
            ## #################################################################
            base <- 31000
            q <- 1
            d <- 700
        } else {
            stop("invalid ipart")
        }
        
        ## parameter settings shared across this set of simulations
        shape <- "line - mean signal in 1st direction"
        k <- 3
        vspikes <- c(1, 10, 50, 100, 1000)
        dfs <- c(-1, 6, 3)
        ni <- 50
        deltas <- seq(2, 10, by=2)
        if (d > 100) { deltas <- 2*deltas }

        ## parse parameter indices by ipart sim index
        ispike <- (floor((ipart - base - 1) / 15) %% 5) + 1
        idelta <- ((ipart - base - 1) %% 5) + 1
        idf <- (floor((ipart - base - 1) / 5) %% 3) + 1

        ## determine simulation parameters
        spike <- vspikes[ispike]
        delta <- deltas[idelta]
        df <- dfs[idf]

        m1 <- c(-delta, rep(0, d-1))
        m2 <- rep(0, d)
        m3 <- c( delta, rep(0, d-1))
        mui <- rbind(m1, m2, m3)        
        vsd <- c(rep(sqrt(spike), q), rep(1, d - q))
    }
    
    ## #########################################################################
    ## K=3 simulation settings [line - mean signal in second direction]
    ## #########################################################################
    if (ipart %in% 32000:32999) {
        if (ipart %in% 32001:32075) {
            ## #################################################################
            ## dim = 1000, spikes = 1
            ## #################################################################
            base <- 32000
            q <- 1
            d <- 700
        } else {
            stop("invalid ipart")
        }
        
        ## parameter settings shared across this set of simulations
        shape <- "line - mean signal in 2nd direction"
        k <- 3
        vspikes <- c(1, 10, 50, 100, 1000)
        dfs <- c(-1, 6, 3)
        ni <- 50
        deltas <- seq(2, 10, by=2)
        if (d > 100) { deltas <- 2*deltas }
        
        ## parse parameter indices by ipart sim index
        ispike <- (floor((ipart - base - 1) / 15) %% 5) + 1
        idelta <- ((ipart - base - 1) %% 5) + 1
        idf <- (floor((ipart - base - 1) / 5) %% 3) + 1

        ## determine simulation parameters
        spike <- vspikes[ispike]
        delta <- deltas[idelta]
        df <- dfs[idf]

        m1 <- c(0, -delta, rep(0, d-2))
        m2 <- rep(0, d)
        m3 <- c(0, delta, rep(0, d-2))
        mui <- rbind(m1, m2, m3)        
        vsd <- c(rep(sqrt(spike), q), rep(1, d - q))
    }
    


    ## #########################################################################
    ## K=3 simulation settings [triangle + jitter - mean signal in first 2 directions]
    ## #########################################################################
    if (ipart %in% 33000:33999) {
        if (ipart %in% 33001:33075) {
            ## #################################################################
            ## dim = 1000, spikes = 1
            ## #################################################################
            base <- 33000
            q <- 1
            d <- 700
        } else {
            stop("invalid ipart")
        }
        
        ## parameter settings shared across this set of simulations
        shape <- "trianglejitter - mean signal in first 2 directions"
        k <- 3
        vspikes <- c(1, 10, 50, 100, 1000)
        dfs <- c(-1, 6, 3)
        ni <- 50
        deltas <- seq(2, 10, by=2)
        if (d > 100) { deltas <- 2*deltas }

        ## parse parameter indices by ipart sim index
        ispike <- (floor((ipart - base - 1) / 15) %% 5) + 1
        idelta <- ((ipart - base - 1) %% 5) + 1
        idf <- (floor((ipart - base - 1) / 5) %% 3) + 1

        ## determine simulation parameters
        spike <- vspikes[ispike]
        delta <- deltas[idelta]
        df <- dfs[idf]

        ## place cluster centers at random jitter from clean "distance" specifications
        m1 <- c(-delta/2, -delta*sqrt(3)/4, rep(0, d-2))
        m2 <- c(      0,  delta*sqrt(3)/4, rep(0, d-2))
        m3 <- c( delta/2, -delta*sqrt(3)/4, rep(0, d-2))
        mui <- rbind(m1, m2, m3) + c(rnorm(k*2, mean=0, sd=delta/4), rep(0, k*(d-2)))
        vsd <- c(rep(sqrt(spike), q), rep(1, d - q))
    }

    ## #########################################################################
    ## K=3 simulation settings [triangle + jitter - mean signal in 2nd and 3rd directions]
    ## #########################################################################
    if (ipart %in% 34000:34999) {
        if (ipart %in% 34001:34075) {
            ## #################################################################
            ## dim = 1000, spikes = 1
            ## #################################################################
            base <- 34000
            q <- 1
            d <- 700
        } else {
            stop("invalid ipart")
        }
        
        ## parameter settings shared across this set of simulations
        shape <- "trianglejitter - mean signal in 2nd and 3rd directions"
        k <- 3
        vspikes <- c(1, 10, 50, 100, 1000)
        dfs <- c(-1, 6, 3)
        ni <- 50
        deltas <- seq(2, 10, by=2)
        if (d > 100) { deltas <- 2*deltas }
        
        ## parse parameter indices by ipart sim index
        ispike <- (floor((ipart - base - 1) / 15) %% 5) + 1
        idelta <- ((ipart - base - 1) %% 5) + 1
        idf <- (floor((ipart - base - 1) / 5) %% 3) + 1
        
        ## determine simulation parameters
        spike <- vspikes[ispike]
        delta <- deltas[idelta]
        df <- dfs[idf]
        
        ## place cluster centers at random jitter from clean "distance" specifications
        m1 <- c(0, -delta/2, -delta*sqrt(3)/4, rep(0, d-3))
        m2 <- c(0,       0,  delta*sqrt(3)/4, rep(0, d-3))
        m3 <- c(0, delta/2, -delta*sqrt(3)/4, rep(0, d-3))
        mui <- rbind(m1, m2, m3) + c(rnorm(k*2, mean=0, sd=delta/4), rep(0, k*(d-2)))
        vsd <- c(rep(sqrt(spike), q), rep(1, d - q))
    }
    
    ## #########################################################################
    ## K=4 simulation settings [tetrahedron + jitter]
    ## #########################################################################
    if (ipart %in% 41000:41999) {
        if (ipart %in% 41001:41075) {
            ## #################################################################
            ## dim = 1000, spikes = 1
            ## #################################################################
            base <- 41000
            q <- 1
            d <- 700
        } else {
            stop("invalid ipart")
        }
        
        ## parameter settings shared across this set of simulations
        shape <- "tetrajitter"
        k <- 4
        vspikes <- c(1, 10, 50, 100, 1000)
        dfs <- c(-1, 6, 3)
        ni <- 50
        deltas <- seq(2, 10, by=2)
        if (d > 100) { deltas <- 2*deltas }

        ## parse parameter indices by ipart sim index
        ispike <- (floor((ipart - base - 1) / 15) %% 5) + 1
        idelta <- ((ipart - base - 1) %% 5) + 1
        idf <- (floor((ipart - base - 1) / 5) %% 3) + 1

        ## determine simulation parameters
        spike <- vspikes[ispike]
        delta <- deltas[idelta]
        df <- dfs[idf]

        ## place cluster centers at random jitter from clean "distance" specifications
        m1 <- c(-delta/2,       0, -delta/sqrt(8), rep(0, d-3))
        m2 <- c( delta/2,       0, -delta/sqrt(8), rep(0, d-3))
        m3 <- c(      0, -delta/2,  delta/sqrt(8), rep(0, d-3))
        m4 <- c(      0,  delta/2,  delta/sqrt(8), rep(0, d-3))
        mui <- rbind(m1, m2, m3, m4) + c(rnorm(k*3, mean=0, sd=delta/4), rep(0, k*(d-3)))
        vsd <- c(rep(sqrt(spike), q), rep(1, d - q))
    }
    
    
    
    ## #########################################################################
    ## K=5 simulation settings [random in K-1 dim sphere]
    ## #########################################################################
    if (ipart %in% 51000:51999) {
        k <- 5
        shape <- "random shell"
        vspikes <- c(1, 10, 50, 100, 1000)
        ni <- 50
        base <- 51000
        q <- 1

        if (ipart %in% 51001:51075) {
            d <- 700
            deltas <- seq(5, 25, by=5)
        } 

        if (((ipart-base-1) %% 15) %in% 0:4) {
            df <- -1
        } else if (((ipart-base-1) %% 15) %in% 5:9) {
            df <- 6
        } else if (((ipart-base-1) %% 15) %in% 10:14) {
            df <- 3
        }
        
        spike <- vspikes[ceiling((ipart - base) / (15))]
        delta <- deltas[((ipart-base-1) %% 5) + 1]
        
        ## place clusters uniform random in k-1 dim sphere w/ radius delta
        mui <- cbind(matrix(rnorm(k*(k-1)), ncol=k, nrow=k-1))
        mui <- scale(mui, center=FALSE) %*% diag(runif(k)^(1/(k-1))) / sqrt(k-2)
        mui <- cbind(t(mui) * delta, matrix(0, nrow=k, ncol=d-k+1))
        vsd <- c(rep(sqrt(spike), q), rep(1, d - q))  
    }
    
    
    
    ## return results
    if (params) {
        output <- list(k=k, ni=ni, d=d, delta=delta, shape=shape,
                       mui=mui, vsd=vsd, df=df)
    } else {
        output <- shc_simulator(k=k, d=d, ni=ni, mui=mui, vsd=vsd, df=df)
    }
    output
}
