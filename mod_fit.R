#==============================================================================#
# Model 1: Iterated Filtering                                                  #
#          Code modified from Elizabeth C. Lee, Andrew S. Azman, and Justin    #
#          Lessler at Johns Hopkins Bloomberg School of Public Health          #
#==============================================================================#

mod_fit <- function(settings, out_dir, period, starts, full_mod) {
  #### make list of parameters being estimated
  if (period == 1) { ## endemic
    est_params <- c(paste0("beta", 1:6), "rho", "tau", "nu")
    rw_sds <- rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02,
                    beta4 = 0.02, beta5 = 0.02, beta6 = 0.02, #sig_sq = 0.02,
                    tau = 0.02, rho = 0.02, nu = 0.02)
    tstart <- 233
    tend <- 430
  } else if (period == 0) { ## epidemic
    est_params <- c(paste0("beta", 1:6), "rho", "tau", "nu", 
                    "E_0", "I_0")
    rw_sds <- rw.sd(beta1 = .02, beta2 = .02, beta3 = .02,
                    beta4 = .02, beta5 = .02, beta6 = .02,
                    tau = 0.02, rho = 0.02, nu = 0.02, #sig_sq = 0.02,
                    E_0 = ivp(0.2), I_0 = ivp(0.2))
    tstart <- 1
    tend <- 232
  } else {
    print("invalid period: must be 0 - epidemic or 1 - endemic")
    return()
  }
  
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
             tic <- Sys.time()
             print(paste(start$parid, "iterated filtering"))
             po <- window(full_mod, start = tstart, end = tend)
             timezero(po) <- tstart - 1
             allpars <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6",
                          "gamma","sigma","theta0","alpha","mu","delta","nu", #"sig_sq",
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
             
             #if (!is.na(ll[1]) & ll[1] > -2000) { ## filter out extreme cases based on log-likelihood
               ## record mean effective sample size
               mean_ess <- vector(length = 10)
               for (i in 1:10) {
                 mean_ess[i] <- mean(pf.lik[[i]]@eff.sample.size)
               }
               
               toc <- Sys.time()
               etime <- toc-tic
               units(etime) <- "hours"
               
               ## simulate
               sims <- sim_data(po, pars = coef(mf.mod), 25)
               ## record simulations
               write_csv(sims, paste0(out_dir, "fitstates_if_", settings$str_if2, 
                                      "_parid", start$parid, "_seed",settings$seed, ".csv"))
               
               ## record parameter estimates
               dummy <- data.frame(model = paste0("if", period), 
                                   parid = start$parid,
                                   as.list(coef(mf.mod)),
                                   loglik = ll[1],
                                   loglik.se = ll[2],
                                   m_ess.min = min(mean_ess),
                                   m_ess.max = max(mean_ess),
                                   etime = as.numeric(etime))
               write_csv(dummy, paste0(out_dir, "parest_if_", settings$str_if2, 
                                       "_parid", start$parid, "_seed", 
                                       settings$seed, ".csv"))
               
               rm(fit_states, fit_states_sims, mf.mod, final_timept)
               gc()
               dummy
             #}
           }
       })  -> prof_if
  return(prof_if)
}
