##
## Given the rds mcmc file calculate the ESS
## Author: Elizabeth Goult  
## Comments:
##

# clear workspace
rm(list = ls())
graphics.off()

# packages
require(FME)
require(mcmcse)
require(tidyverse)


# function
CalcMESS<- function(mcmc){
    
    n_par<-length(mcmc$bestpar)
    # calculate minimum ESS for problem at 5%    
    min_val <- minESS(p = n_par)

    # calculate the ess
    m_ess<-multiESS(x = mcmc$pars)

    min_eps <- minESS(p = n_par,  ess = m_ess)

    return(c(miness = min_val, ess = m_ess, epsilon = min_eps, enough = m_ess > min_val, npar = n_par, niter = nrow(mcmc$pars)))
}

# read in mcmc rds object
mcmc.path<- "results/results_20211030/petar/"
mcmc.strings<- paste0(mcmc.path, list.files(path = mcmc.path) )
mcmc <- lapply(mcmc.strings, readRDS)


MESS <- lapply(mcmc, CalcMESS) %>%
    bind_rows() %>%
    bind_cols(run = mcmc.strings)


write.csv(MESS, paste0(mcmc.path, "effective_sample_sizes.csv"))
