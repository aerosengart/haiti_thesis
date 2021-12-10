# Loading the necessary libraries
library(tidyverse)
library(doParallel)
library(foreach)
library(panelPomp)
library(doRNG)
library(tictoc)

cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
# cores <- 20
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)

# Set Run Level for debug, timing, and full computation
RUN_LEVEL = 3

chol_Np <-           switch(RUN_LEVEL, 50, 1000, 5000)  # Number of particle filters
chol_Nmif <-         switch(RUN_LEVEL,  3, 15, 100)  # Number of MIF iterations
# chol_Nreps_eval <-   switch(RUN_LEVEL,  2,  10,  20)  # Number of times pfilter will be run to estimate likelihood
chol_Nreps_local <-  switch(RUN_LEVEL,  2,  10, 36)  # Number of times to run MIF at "MLE"

h1_panel <- haiti1_panel()

rws <- rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02,
             beta4 = 0.02, beta5 = 0.02, beta6 = 0.02,
             rho_epi = ifelse(time > 232, 0.0, 0.02),
             rho_end = ifelse(time <= 232, 0.0, 0.02),
             tau_epi = ifelse(time > 232, 0.0, 0.02),
             tau_end = ifelse(time <= 232, 0.0, 0.02),
             sig_sq_epi = ifelse(time > 232, 0.0, 0.02),
             sig_sq_end = ifelse(time <= 232, 0.0, 0.02),
             nu = 0.02,
             E_0 = ifelse(time > 232, 0.0, ivp(0.2)),
             I_0 = ifelse(time > 232, 0.0, ivp(0.2)))

# Get all of the cores available.
registerDoRNG(1851563)

# Local MIF at "MLE"
stew(file = 'model1/panel/mif_panel.rda', {

  t2 <- system.time({
    # Run local MIF chol_Nreps_local times
    foreach(
      i=1:chol_Nreps_local,
      .packages = c('panelPomp'),
      .combine = c
    ) %dopar% {
      mif2(
        h1_panel,
        Np = chol_Np,
        Nmif = chol_Nmif,
        cooling.fraction.50 = 0.5,
        rw.sd = rws,
        cooling.type = 'hyperbolic'
      )
    } -> m2
  })
})


