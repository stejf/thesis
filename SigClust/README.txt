= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
Statistical Significance of Clustering (SigClust) simulations code
Author: Steve Fayek (steve.fayek20@imperial.ac.uk)
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

This directory contains two folders:
    * multiple_starts_experiment : code used for experiments to justify the use of 
				   just one random start for sigclust 
    * simulations		 : contains the code for the true, sample, hard 
				   and soft methods, for three distributions:
				   single Gaussian and mixture of two Gaussians 
				   with mean signal in one direction and in all 
				   directions.
				   
The code underlying the sigclust() function can be accessed at:
https://rdrr.io/cran/sigclust/src/R/sigclust.R
