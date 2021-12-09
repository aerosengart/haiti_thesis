## load libraries
library(plyr)
library(tidyverse)
library(pomp)
library(devtools)
library(data.table)
library(foreach)
library(doParallel)
library(doRNG)
library(iterators)
library(pkgbuild)
library(gridExtra)
library(GGally)
library(subplex)
library(dplyr)
cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)

##----settings --- change if desired--------------------------------------------
# Set Run Level
RUN_LEVEL = 1
num_parts <- switch(RUN_LEVEL, 50, 1e3, 5e3)  # Number of particles
num_iters <- switch(RUN_LEVEL, 5, 100, 200)  # Number of MIF iterations
num_profs <- switch(RUN_LEVEL,  1,  3,  10)  # Number of profiles
num_sims <- switch(RUN_LEVEL,  25,  50,  100)  # Number of times pfilter will be run to estimate likelihood

####EPIDEMIC MODEL FITTING######################################################
#### settings and data
fix_epi_settings <- function(nprofs = 3, nmif = 100, nparticles = 1000, seed = 20190415) {
  epi_nweeks <- 232 ## 180 end of 3/2014; 219 end of 2014; 232 end of 3/2015
  epi_seed <- seed ## seed
  epi_nprofs <- nprofs ## number of profiles in global search
  epi_nmif <- nmif ## number of iterations of mif2
  epi_nparticles <- nparticles ## number of particles for mif2 and pfilter
  modcode <- "final"
  
  ## create strings for filenames
  epi_str_if <- paste0("epi", epi_nweeks, paste0("_nprofs", epi_nprofs),
                       "_nmif", epi_nmif, "_nparticles", epi_nparticles,
                       "_seed", epi_seed)
  epi_str_if2 <- paste0("epi", epi_nweeks, paste0("_nprofs", epi_nprofs),
                        "_nmif", epi_nmif, "_nparticles", epi_nparticles)
  
  return(list(nweeks = epi_nweeks, seed = epi_seed, nmif = epi_nmif,
              nparticles = epi_nparticles, mcode = modcode,
              nprofs = epi_nprofs, str_if = epi_str_if, str_if2 = epi_str_if2))
}
episettings <- fix_epi_settings(nprofs = num_profs, nparticles = num_parts, nmif = num_iters)

fix_end_settings <- function(episettings, seed = 4172019) {
  end_seed <- 4172019
  end_nmif <- episettings$nmif
  end_nparticles <- episettings$nparticles
  modcode <- paste0(episettings$mcode, "1")
  
  end_str_if <- paste0("end", episettings$nweeks, "_nmif", end_nmif,
                       "_nparticles", end_nparticles, "_seed", end_seed)
  end_str_if2 <- paste0("end", episettings$nweeks, "_nmif", end_nmif,
                        "_nparticles", end_nparticles)
  
  return(list(seed=end_seed, nmif=end_nmif, nparticles=end_nparticles,
              mcode=modcode, str_if = end_str_if, str_if2 = end_str_if2))
}
endsettings <- fix_end_settings(episettings)

generate_epi_params <- function(out_dir) {
  parseed <- fix_epi_settings()$seed
  nweeks <- fix_epi_settings()$nweeks
  mcode <- fix_epi_settings()$mcode
  nprofs <- fix_epi_settings()$nprofs
  
  if(file.exists(paste0(out_dir, "starts_", mcode, "_epi", nweeks,
                        paste0("_nprofs", nprofs), "_", parseed, ".csv"))) {
    starts_epi <- read_csv(paste0(out_dir, "starts_", mcode, "_epi", nweeks,
                                  paste0("_nprofs", nprofs), "_", parseed, ".csv"))
  }
  else{
    set.seed(parseed)
    pop0 <- get_haiti_pop()
    
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
    E0 <- 10/pop0 ## rpois(nsamps, 10)/pop
    I0 <- 10/pop0 ## rpois(nsamps, 10)/pop
    A0 <- 0.0/pop0
    R0 <- 0.000
    S0 <- 1-R0-I0-E0-A0
    
    ## beta parameter settings
    blo <- 1E-9; bup <- 10 ## uniform beta settings
    ## median R0 among init values: median(exp(rnorm(1000, log(bmn), bse)))*2/7
    bmn <- 4.5; bse <- 0.5
    ## nu parameter settings
    nlo <- 0.95; nup <- 1
    
    rhosamps <- profile_design(rho = seq(1E-8, 1, length = 30),
                               upper = c(tau = 20, beta1 = bup, beta2 = bup,
                                         beta3 = bup, beta4 = bup, beta5 = bup,
                                         beta6 = bup, nu = nup),
                               lower = c(tau = 1, beta1 = blo, beta2 = blo,
                                         beta3 = blo, beta4 = blo, beta5 = blo,
                                         beta6 = blo, nu = nlo),
                               nprof = nprofs)
    tausamps <- profile_design(tau = seq(1, 20, length = 30),
                               upper = c(rho = 1, beta1 = bup, beta2 = bup,
                                         beta3 = bup, beta4 = bup, beta5 = bup,
                                         beta6 = bup, nu = nup),
                               lower = c(rho = 1E-8, beta1 = blo, beta2 = blo,
                                         beta3 = blo, beta4 = blo, beta5 = blo,
                                         beta6 = blo, nu = nlo),
                               nprof = nprofs)
    betasamps <- profile_design(beta1 = seq(blo, bup, length = 30),
                                upper = c(rho = 1, tau = 20, beta2 = bup,
                                          beta3 = bup, beta4 = bup, beta5 = bup,
                                          beta6=bup, nu = nup),
                                lower = c(rho = 1E-8, tau = 1, beta2 = blo,
                                          beta3=blo, beta4 = blo, beta5 = blo,
                                          beta6=blo, nu = nlo),
                                nprof = nprofs)
    nusamps <- profile_design(nu = seq(nlo, nup, length = 10),
                              upper = c(rho = 1, tau = 20, beta1 = bup, beta2 = bup,
                                        beta3 = bup, beta4 = bup, beta5 = bup,
                                        beta6 = bup),
                              lower = c(rho = 1E-8, tau = 1, beta1 = blo,
                                        beta2 = blo, beta3 = blo, beta4 = blo,
                                        beta5 = blo, beta6 = blo),
                              nprof = nprofs)
    starts_epi <- bind_rows(rhosamps, tausamps, betasamps, nusamps) %>%
      mutate(parid = seq_along(rho)) %>%
      mutate(gamma = gamma, sigma = sigma, theta0 = theta0, alpha = alpha,
             mu = mu, delta = delta, nu = nu, S_0 = S0, E_0 = E0, I_0 = I0,
             A_0 = A0, R_0 = R0, pop_0 = pop0) %>%
      select(parid, rho, tau, beta1, beta2, beta3, beta4, beta5,
             beta6, nu, gamma, sigma, theta0, alpha, mu, delta,
             S_0, E_0, I_0, A_0, R_0, pop_0)
  }
  return(starts_epi)
}
starts <- generate_epi_params(out_dir = out_dir) ## generate starting values for parameter estimation

pop.haiti <- sum(c(746236, 1727524, 4029705, 728807, 1067177, 774976, 342525, 393967, 632601, 468301)) ## get haiti population (country-wide)

get_mspp_agg_data <- function() {
  ## Saturday dates represent last day of week cumulatively
  allDat <- read_csv("Data/haiti-data-from-2010-10-to-2019-01.csv", skip = 1,
                     col_names = c("date_sat_orig", "report", "Artibonite",
                                   "Centre", "Grand_Anse", "Nippes", "Nord",
                                   "Nord_Est", "Nord_Ouest", "Ouest", "Sud",
                                   "Sud_Est"), col_types = "cciiiiiiiiii")
  
  splitDate <- strsplit(allDat$date_sat_orig, "-")
  setattr(splitDate[[1]], 'names', c("year", "month", "day"))
  dateDf <- tibble::as_tibble(as.data.frame(do.call(rbind, splitDate))) %>%
    mutate(month = as.character(month)) %>%
    mutate(day = as.character(day)) %>%
    mutate(year = as.character(year)) %>%
    mutate(month = ifelse(nchar(month) == 1, paste0("0", month), month)) %>%
    mutate(day = ifelse(nchar(day) == 1, paste0("0", day), day)) %>%
    mutate(date_sat = as.Date(paste(year, month, day, sep = "-"), origin = "1900-01-01"))
  
  fullDateVec <- data.frame(date_sat = seq(min(dateDf$date_sat), max(dateDf$date_sat), by = 7))
  
  cleanDat <- allDat %>%
    mutate(date_sat = as.Date(dateDf$date_sat, origin = "1900-01-01")) %>%
    select(-date_sat_orig) %>%
    full_join(fullDateVec, by = c("date_sat")) %>%
    arrange(date_sat) %>%
    mutate(week = seq_along(date_sat)) %>%
    tidyr::gather(department, cases, Artibonite:Sud_Est)
  
  aggDat <- cleanDat %>%
    group_by(week) %>% ## "day" or "week"
    summarise(cases = sum(cases, na.rm=TRUE))
  
  return(aggDat)
}
haiti.dat <- get_mspp_agg_data() ## get case data (country-wide)

haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= episettings$nweeks) ## get case data (epidemic period)
haiti.dat.end <- haiti.dat %>% dplyr::filter(week > episettings$nweeks) %>% ## get case data (endemic period)
  dplyr::mutate(week_end = seq_along(week))

make_covartab <- function(t0, tmax, byt = 1, nbasis = 6, degree = 6, per = 52.14) {
  tbasis <- seq(from = t0, to = tmax, by = byt)
  covartab <-  data.frame(cbind(time = tbasis,
                                periodic.bspline.basis(x = tbasis,
                                                       nbasis = nbasis,
                                                       degree = degree,
                                                       period = per,
                                                       names = "seas%d")))
  
  return(covartab)
}
covartab <- make_covartab(0, nrow(haiti.dat)+52*11+1, byt=1, degree=6, nbasis=6, per=52.14)

source("mod_build.R")
haiti.mod <- mod_build(dat = haiti.dat, covar = covartab) ## build model

#### perform mif2 on epidemic period
mod_fit <- function(settings, out_dir, period, starts, full_mod) {
  #### make list of parameters being estimated
  if (period == 1) { ## endemic
    est_params <- c(paste0("beta", 1:6), "rho", "tau", "nu")
    rw_sds <- rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02,
                    beta4 = 0.02, beta5 = 0.02, beta6 = 0.02,
                    tau = 0.02, rho = 0.02, nu = 0.02)
    tstart <- 233
    tend <- 430
  } else if (period == 0) { ## epidemic
    est_params <- c(paste0("beta", 1:6), "rho", "tau", "nu", 
                    "E_0", "I_0")
    rw_sds <- rw.sd(beta1 = .02, beta2 = .02, beta3 = .02,
                    beta4 = .02, beta5 = .02, beta6 = .02,
                    tau = 0.02, rho = 0.02, nu = 0.02,
                    E_0 = ivp(0.2), I_0 = ivp(0.2))
    tstart <- 1
    tend <- 232
  } else {
    print("invalid period: must be 0 - epidemic or 1 - endemic")
    return()
  }
  
  bake(file = "epi_if_script_output.rda"), 
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
                          "gamma","sigma","theta0","alpha","mu","delta","nu",
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
             
             ## get likelihood estimate 
             pf.lik <- replicate(10, pfilter(mf.mod, Np = settings$nparticles)) 
             ll <- sapply(pf.lik, logLik)
             ll <- logmeanexp(ll, se = TRUE)
             
             if (!is.na(ll[1]) & ll[1] > -2000) { ## filter out extreme cases based on log-likelihood
               ## record mean effective sample size
               mean_ess <- vector(length = 10)
               for (i in 1:10) {
                 mean_ess[i] <- mean(pf.lik[[i]]@eff.sample.size)
               }
               
               toc <- Sys.time()
               etime <- toc-tic
               units(etime) <- "hours"
               
               ## record parameter estimates
               dummy <- data.frame(model = paste0("if", period), 
                                   parid = start$parid,
                                   as.list(coef(mf.mod)),
                                   loglik = ll[1],
                                   loglik.se = ll[2],
                                   m_ess.min = min(mean_ess),
                                   m_ess.max = max(mean_ess),
                                   etime = as.numeric(etime))
               
               rm(fit_states, fit_states_sims, mf.mod, final_timept)
               gc()
               dummy
             }
           }
       }) -> prof_if
  return(prof_if)
}
prof_if <- mod_fit(episettings, out_dir, 0, starts, haiti.mod)
