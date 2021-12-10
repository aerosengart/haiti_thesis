library(pomp)
library(haitipkg)
library(foreach)
library(doParallel)
library(doRNG)
library(tidyverse)
library(dplyr)
library(data.table)

source("model1/mod_tools.R")
source("model1/mod_fit.R")

cores <- as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset = NA))
if(is.na(cores)) cores <- 36
registerDoParallel(cores)

## run level for testing/computing
RUN_LEVEL = 3

# Set Run Level
num_parts <- switch(RUN_LEVEL, 50, 1e3, 5e3)  # Number of particles
num_iters <- switch(RUN_LEVEL, 5, 100, 200)  # Number of MIF iterations
num_profs <- switch(RUN_LEVEL,  1,  3,  5)  # Number of profiles

## make model
full_mod <- haiti1()
# base_dir <- "model1/adj_mod/"
# out_dir <- "model1/adj_mod/end/"
# read_dir <- "model1/adj_mod/epi/"
base_dir <- "report_data/adj_mod/"
out_dir <- "report_data/adj_mod/end/"
read_dir <- "report_data/adj_mod/epi/"
dir.create(base_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)
episettings <- fix_epi_settings(nprofs = num_profs, nparticles = num_parts, nmif = num_iters)
endsettings <- fix_end_settings(episettings)
full_mod <- full_mod %>%
  window(start = episettings$nweeks + 1, end = 430)

pop_haiti <- get_haiti_pop()
haiti.dat <- haiti1_agg_data() ## get case data (country-wide)

set.seed(endsettings$seed)

starts <- generate_end_params(out_dir_epi = read_dir, out_dir_end = out_dir,
                              episettings = episettings, endsettings = endsettings)

write_csv(starts, paste0(out_dir, "end_starts.csv"))

settings <- endsettings
est_params <- c(paste0("beta", 1:6), "rho", "tau", "nu", "sig_sq")
rw_sds <- rw.sd(beta1 = .02, beta2 = .02, beta3 = .02,
                beta4 = .02, beta5 = .02, beta6 = .02,
                tau = 0.02, rho = 0.02, nu = 0.02,
                sig_sq = 0.02)

bake(file = paste0(out_dir, "/iffits_", settings$str_if, ".rds"),
   seed = settings$seed, {
     starts_if <- starts
     foreach(start = iter(starts_if, by = 'row'),
             .combine = rbind, .inorder = FALSE,
             .packages = c("pomp", "magrittr"),
             .errorhandling = c('remove'),
             .export = c("haiti.dat", "num_betas", "pop.haiti",
                         "covartab", "settings"),
             .noexport = c(),
             .verbose = TRUE) %dopar%
       {
         po <- full_mod
         timezero(po) <- po@times[1] - 1
         allpars <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6",
                      "gamma","sigma","theta0","alpha","mu","delta","nu", "sig_sq",
                      "S_0","E_0","I_0","A_0","R_0", "pop_0")
         coef(po) <- unlist(start[which(names(start) %in% allpars)])

         ## if no Exposed/Infectious, make it small but nonzero
         if (coef(po, "E_0") == 0.0 ) {
           coef(po, "E_0") <- 1e-9
         }
         if (coef(po, "I_0") == 0.0) {
           coef(po, "I_0") <- 1e-9
         }

         ## perform iterated filtering
         mf.mod <- mif2(po, Nmif = settings$nmif,
                        rw.sd = rw_sds,
                        Np = settings$nparticles,
                        cooling.type = "hyperbolic",
                        cooling.fraction.50 = 0.5,
                        verbose = FALSE)
         save(mf.mod, file = paste0(out_dir, "fitted_mod_", start$parid, ".rda"))

         ## get likelihood estimate
         pf.lik <- replicate(10, pfilter(mf.mod, Np = settings$nparticles))
         ll <- sapply(pf.lik, logLik)
         ll <- logmeanexp(ll, se = TRUE)

         #if (!is.na(ll[1])) { ## filter out extreme cases based on log-likelihood
           ## record mean effective sample size
           mean_ess <- vector(length = 10)
           for (i in 1:10) {
             mean_ess[i] <- mean(pf.lik[[i]]@eff.sample.size)
           }

           ## simulate
           sims <- sim_data(po, pars = coef(mf.mod), 25)
           ## record simulations
           write_csv(sims, paste0(out_dir, "fitstates_if_parid", start$parid, ".csv"))

           ## record parameter estimates
           dummy <- data.frame(model = paste0("if"),
                               parid = start$parid,
                               as.list(coef(mf.mod)),
                               loglik = ll[1],
                               loglik.se = ll[2],
                               m_ess.min = min(mean_ess),
                               m_ess.max = max(mean_ess))
           write_csv(dummy, paste0(out_dir, "parest_if_", settings$str_if2,
                                   "_parid", start$parid, "_seed",
                                   settings$seed, ".csv"))

           rm(mf.mod, pf.lik, sims)
           gc()
           dummy
          #}
   }
})  -> prof_if

save(prof_if, file = "model1/adj_mod/end/adj_end_mif.rda")
