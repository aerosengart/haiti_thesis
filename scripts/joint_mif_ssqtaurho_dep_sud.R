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
num_profs <- switch(RUN_LEVEL,  1,  3,  3)  # Number of profiles

## make model
full_mod <- haiti1_dep("Sud", "id0")
out_dir <- "model1/joint_mif/str_sud"
dir.create(out_dir, showWarnings = FALSE)
episettings <- fix_epi_settings(nprofs = num_profs, nparticles = num_parts, nmif = num_iters)

#pop_haiti <- get_haiti_pop()
pop_haiti <- 632601
haiti.dat <- haiti1_data() ## get case data (country-wide)

set.seed(episettings$seed)

## Fixed parameters
gamma <- 7/2 ## 2 day infectious period
sigma <- 7/1.4 ## 1.4 day latent period
theta0 <- 0 ## 0% going to asymptomatic in base model
alpha <- 7/2920 ## 8 year mean duration of natural immunity
## 22.6/1000 average annual birth rate, adjusted for compounding by week
mu <- ((1+22.6/1000)^(1/52.14))-1
## 7.5/1000 average annual death rate, adjusted for compounding by week
delta <- ((1+7.5/1000)^(1/52.14))-1

## starting states
E0 <- 10/pop_haiti ## rpois(nsamps, 10)/pop
I0 <- 10/pop_haiti ## rpois(nsamps, 10)/pop
A0 <- 0.0
R0 <- 0.0
S0 <- 1-R0-I0-E0-A0

## beta parameter settings
blo <- 1E-9; bup <- 10 ## uniform beta settings
## median R0 among init values: median(exp(rnorm(1000, log(bmn), bse)))*2/7
bmn <- 4.5; bse <- 0.5
## nu parameter settings
nlo <- 0.95; nup <- 1
## sigma_se parameter settings
sigsqlo <- 1E-9; sigsqup <- 5

rhosamps <- profile_design(rho = seq(1E-8, 1, length = 30),
                           upper = c(tau = 20, beta1 = bup, beta2 = bup,
                                     beta3 = bup, beta4 = bup, beta5 = bup,
                                     beta6 = bup, nu = nup, sig_sq = sigsqup),
                           lower = c(tau = 1, beta1 = blo, beta2 = blo,
                                     beta3 = blo, beta4 = blo, beta5 = blo,
                                     beta6 = blo, nu = nlo, sig_sq = sigsqlo),
                           nprof = num_profs)
tausamps <- profile_design(tau = seq(1, 20, length = 30),
                           upper = c(rho = 1, beta1 = bup, beta2 = bup,
                                     beta3 = bup, beta4 = bup, beta5 = bup,
                                     beta6 = bup, nu = nup, sig_sq = sigsqup),
                           lower = c(rho = 1E-8, beta1 = blo, beta2 = blo,
                                     beta3 = blo, beta4 = blo, beta5 = blo,
                                     beta6 = blo, nu = nlo, sig_sq = sigsqlo),
                           nprof = num_profs)
betasamps <- profile_design(beta1 = seq(blo, bup, length = 30),
                            upper = c(rho = 1, tau = 20, beta2 = bup,
                                      beta3 = bup, beta4 = bup, beta5 = bup,
                                      beta6=bup, nu = nup, sig_sq = sigsqup),
                            lower = c(rho = 1E-8, tau = 1, beta2 = blo,
                                      beta3=blo, beta4 = blo, beta5 = blo,
                                      beta6=blo, nu = nlo, sig_sq = sigsqlo),
                            nprof = num_profs)
nusamps <- profile_design(nu = seq(nlo, nup, length = 30),
                          upper = c(rho = 1, tau = 20, beta1 = bup, beta2 = bup,
                                    beta3 = bup, beta4 = bup, beta5 = bup,
                                    beta6 = bup, sig_sq = sigsqup),
                          lower = c(rho = 1E-8, tau = 1, beta1 = blo,
                                    beta2 = blo, beta3 = blo, beta4 = blo,
                                    beta5 = blo, beta6 = blo, sig_sq = sigsqlo),
                          nprof = num_profs)
sigsqsamps <- profile_design(sig_sq = seq(sigsqlo, sigsqup, length = 30),
                             upper = c(rho = 1, tau = 20, beta1 = bup, beta2 = bup,
                                       beta3 = bup, beta4 = bup, beta5 = bup,
                                       beta6 = bup, nu = nup),
                             lower = c(rho = 1E-8, tau = 1, beta1 = blo,
                                       beta2 = blo, beta3 = blo, beta4 = blo,
                                       beta5 = blo, beta6 = blo, nu = nlo),
                             nprof = num_profs)
starts <- bind_rows(rhosamps, tausamps, betasamps, nusamps, sigsqsamps) %>%
  mutate(parid = seq_along(rho)) %>%
  mutate(theta0 = theta0, mu = mu, delta = delta, nu = nu, sigma = sigma, alpha = alpha, kappa = 0,
         gamma = gamma, S_0 = S0, E_0 = E0, I_0 = I0, A_0 = A0, R_0 = R0, pop_0 = pop_haiti) %>%
  select(parid, rho, tau, beta1, beta2, beta3, beta4, beta5,
         beta6, nu, gamma, sigma, theta0, alpha, mu, delta, sig_sq, kappa,
         S_0, E_0, I_0, A_0, R_0, pop_0)
add_starts <- starts %>%
  dplyr::select(rho, tau, sig_sq) %>%
  dplyr::rename(rho_end = rho, sig_sq_end = sig_sq, tau_end = tau)
starts <- cbind(starts, add_starts) %>%
  dplyr::rename(rho_epi = rho, tau_epi = tau, sig_sq_epi = sig_sq)

write_csv(starts, paste0(out_dir, "/starts", paste0("_nprofs", num_profs), "_",
                             episettings$parseed, ".csv"))

settings <- episettings
est_params <- c(paste0("beta", 1:6),
                "rho_epi", "rho_end", "tau_epi",
                "tau_end", "nu", "sig_sq_epi", "sig_sq_end", "E_0", "I_0")
rw_sds <- rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02,
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

bake(file = paste0(out_dir, "/iffits_joint.rds"),
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
         timezero(po) <- 0
         allpars <- c("rho_epi", "rho_end",
                      "beta1", "beta2", "beta3",
                      "beta4", "beta5", "beta6",
                      "sig_sq_epi", "sig_sq_end", "tau_epi", "tau_end",
                      "gamma","sigma","theta0","alpha","mu","delta","nu", "kappa",
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
         mf.mod <- mif2(po,
                        #Nmif = 50,
                        Nmif = settings$nmif,
                        rw.sd = rw_sds,
                        #Np = 1000,
                        Np = settings$nparticles,
                        cooling.type = "hyperbolic",
                        cooling.fraction.50 = 0.5,
                        verbose = FALSE)

         epi_mod <- window(mf.mod, start = 1, end = 232)
         timezero(epi_mod) <- 0
         end_mod <- window(mf.mod, start = 233, end = 430)
         timezero(end_mod) <- 232

         ## get likelihood estimate
         epi.lik <- replicate(10, pfilter(epi_mod, Np = settings$nparticles))
         end.lik <- replicate(10, pfilter(end_mod, Np = settings$nparticles))
         full.lik <- replicate(10, pfilter(mf.mod, Np = settings$nparticles))
         epi_ll <- sapply(epi.lik, logLik)
         epi_ll <- logmeanexp(epi_ll, se = TRUE)
         end_ll <- sapply(end.lik, logLik)
         end_ll <- logmeanexp(end_ll, se = TRUE)
         full_ll <- sapply(full.lik, logLik)
         full_ll <- logmeanexp(full_ll, se = TRUE)

	save(mf.mod, file = paste0(out_dir, "/fitted_mod_", start$parid, ".rda"))

        # if (!is.na(epi_ll[1]) && !is.na(end_ll[1])) { ## filter out extreme cases based on log-likelihood
           ## simulate
           sims <- sim_data(po, pars = coef(mf.mod), 25)
           ## record simulations
           write_csv(sims, paste0(out_dir, "/fitstates_if_parid", start$parid, ".csv"))

           ## record parameter estimates
           dummy <- data.frame(model = paste0("if"),
                               parid = start$parid,
                               as.list(coef(mf.mod)),
                               epi_loglik = epi_ll[1],
                               epi_loglik.se = epi_ll[2],
                               end_loglik = end_ll[1],
                               end_loglik.se = end_ll[2],
                               full_loglik = full_ll[1],
                               full_loglik.se = full_ll[2])
           write_csv(dummy, file = paste0(out_dir, "/parest_if_parid", start$parid, ".csv"))

           rm(mf.mod, epi_mod, end_mod, epi.lik, end.lik, sims)
           gc()
           dummy
         # }
   }
})  -> prof_if

save(prof_if, file = "model1/joint_mif/str_sud/str_sud.rda")
