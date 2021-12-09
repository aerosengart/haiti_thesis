#==============================================================================#
# Model 1: Trajectory Matching                                                 #
#          Code modified from Elizabeth C. Lee, Andrew S. Azman, and Justin    #
#          Lessler at Johns Hopkins Bloomberg School of Public Health          #
#==============================================================================#

mod_traj <- function(settings, out_dir, period, starts, full_mod) {
  #### make list of parameters being estimated
  if (period == 1) { ## endemic
    est_params <- c(paste0("beta", 1:6), "rho", "tau") ## removed NU
    tstart <- 233
    tend <- 430
  } else if (period == 0) { ## epidemic
    est_params <- c(paste0("beta", 1:6), "rho", "tau",
                    "E_0", "I_0")
    tstart <- 1
    tend <- 232
  } else {
    print("invalid period: must be 0 - epidemic or 1 - endemic")
    return()
  }

  ## -------------------------------------- ##
  ## perform epi trajectory matching fits
  bake(file = paste0(out_dir, "/tmfits.rds"),
       seed = settings$seed, {
    foreach(start = iter(starts, by = 'row'),
            .combine = rbind, .inorder = FALSE,
            .packages = c("pomp", "magrittr"),
            .errorhandling = c('remove'),
            .export = c("haiti.dat", "pop.haiti", "covartab", "settings"),
            .noexport = c(),
            .verbose = TRUE) %dopar%
      {
        tm <- window(full_mod, start = tstart, end = tend)
        allpars <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6",
                     "gamma","sigma","theta0","alpha","mu","delta","nu",
                     "S_0","E_0","I_0","A_0","R_0", "pop_0")
        coef(tm) <- unlist(start[which(names(start) %in% allpars)])

        # create trajectory matching objective function
        traj_func <- traj_objfun(tm,
                                 est = est_params,
                                 params = coef(tm),
                                 dmeasure = tm@dmeasure,
                                 partrans = tm@partrans)
        est_params_init <- tm@params
        est_params_init <- est_params_init[which(names(est_params_init) %in% est_params)]

        # do trajectory matching
        traj_fit <- optim(par = est_params_init, fn = traj_func)
        traj_func(traj_fit$par)

        # estimates to natural scale
        est_params <- coef(traj_func)

        # likelihood of estimates
        loglik_est <- logLik(traj_func)

        dummy <- data.frame(model = "tm",
                            parid = start$parid,
                            as.list(est_params),
                            loglik = loglik_est,
                            conv = traj_fit$convergence)
        write_csv(dummy, paste0(out_dir, "parest_tm_parid", start$parid, "_seed", settings$seed, ".csv"))

        rm(tm, traj_fit)
        gc()
        dummy
      }
  })  -> prof_tm
  return(prof_tm)
}



