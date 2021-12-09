#==============================================================================#
# Model 1: Helper functions for data, model building, and model assessment     #
#          Code modified from Elizabeth C. Lee, Andrew S. Azman, and Justin    #
#          Lessler at Johns Hopkins Bloomberg School of Public Health          #
#==============================================================================#

##############################################################################--
####  DATA/SET-UP FUNCTIONS                                                 ####
##############################################################################--
#### get population of Haiti at t0 (use sum of department populations)
get_haiti_pop <- function(departement = 'Total') {
  ## population by departement
  pops <- c(746236, 1727524, 4029705, 728807, 1067177,
            774976, 342525, 393967, 632601, 468301) %>%
    as.data.frame()
  rownames(pops) <- c("Artibonite", "Centre", "Grande_Anse", "Nippes",
                      "Nord", "Nord_Est", "Nord_Ouest", "Ouest", "Sud", "Sud_Est")
  if (departement == 'Total') {
    pop <- sum(pops)
  } else {
    pop <- pops[departement, ]
  }
  return(pop)
}

#### load mspp data aggregated to country scale
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

#### load mspp data aggregated to department scale
get_mspp_dept_data <- function() {
  ## Saturday dates represent last day of week cumulatively
  allDat <- read_csv("Data/haiti-data-from-2010-10-to-2019-01.csv", skip = 1,
                     col_names = c("date_sat_orig", "report", "Artibonite",
                                   "Centre", "Grand_Anse", "Nippes", "Nord",
                                   "Nord_Est", "Nord_Ouest", "Ouest", "Sud",
                                   "Sud_Est"), col_types = "cciiiiiiiiii")

  splitDate <- strsplit(allDat$date_sat_orig, "-")
  setattr(splitDate[[1]], 'names', c("year", "month", "day"))
  dateDf <- as_tibble(as.data.frame(do.call(rbind, splitDate))) %>%
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

  return(cleanDat)
}

#### construct seasonal b-spline terms
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

##############################################################################--
####  UTILITY FUNCTIONS FOR MODEL BUILDING                                  ####
##############################################################################--
#### settings for epidemic model
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

#### settings for the endemic model
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

#### set parameters for epidemic fit model
generate_epi_params <- function(out_dir, departement = 'Total') {
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
    pop0 <- get_haiti_pop(departement)

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
    siglo <- 1E-9; sigup <- 20

    rhosamps <- profile_design(rho = seq(1E-8, 1, length = 15),
                               upper = c(tau = 20, beta1 = bup, beta2 = bup,
                                         beta3 = bup, beta4 = bup, beta5 = bup,
                                         beta6 = bup, nu = nup, sig_sq = sigup),
                               lower = c(tau = 1, beta1 = blo, beta2 = blo,
                                         beta3 = blo, beta4 = blo, beta5 = blo,
                                         beta6 = blo, nu = nlo, sig_sq = siglo),
                               nprof = nprofs)
    tausamps <- profile_design(tau = seq(1, 20, length = 15),
                               upper = c(rho = 1, beta1 = bup, beta2 = bup,
                                         beta3 = bup, beta4 = bup, beta5 = bup,
                                         beta6 = bup, nu = nup, sig_sq = sigup),
                               lower = c(rho = 1E-8, beta1 = blo, beta2 = blo,
                                         beta3 = blo, beta4 = blo, beta5 = blo,
                                         beta6 = blo, nu = nlo, sig_sq = siglo),
                               nprof = nprofs)
    betasamps <- profile_design(beta1 = seq(blo, bup, length = 15),
                                upper = c(rho = 1, tau = 20, beta2 = bup,
                                          beta3 = bup, beta4 = bup, beta5 = bup,
                                          beta6=bup, nu = nup, sig_sq = sigup),
                                lower = c(rho = 1E-8, tau = 1, beta2 = blo,
                                          beta3=blo, beta4 = blo, beta5 = blo,
                                          beta6=blo, nu = nlo, sig_sq = siglo),
                                nprof = nprofs)
    nusamps <- profile_design(nu = seq(nlo, nup, length = 15),
                              upper = c(rho = 1, tau = 20, beta1 = bup, beta2 = bup,
                                        beta3 = bup, beta4 = bup, beta5 = bup,
                                        beta6 = bup, sig_sq = sigup),
                              lower = c(rho = 1E-8, tau = 1, beta1 = blo,
                                        beta2 = blo, beta3 = blo, beta4 = blo,
                                        beta5 = blo, beta6 = blo, sig_sq = siglo),
                              nprof = nprofs)
    sigsamps <- profile_design(sig_sq = seq(siglo, sigup, length = 15),
                              upper = c(rho = 1, tau = 20, beta1 = bup, beta2 = bup,
                                        beta3 = bup, beta4 = bup, beta5 = bup,
                                        beta6 = bup, nu = nup),
                              lower = c(rho = 1E-8, tau = 1, beta1 = blo,
                                        beta2 = blo, beta3 = blo, beta4 = blo,
                                        beta5 = blo, beta6 = blo, nu = nlo),
                              nprof = nprofs)
    starts_epi <- bind_rows(rhosamps, tausamps, betasamps, nusamps, sigsamps) %>%
      mutate(parid = seq_along(rho)) %>%
      mutate(theta0 = theta0, alpha = alpha,
             mu = mu, delta = delta, nu = nu, gamma = gamma, sigma = sigma,
             sig_sq = sig_sq, S_0 = S0, E_0 = E0, I_0 = I0,
             A_0 = A0, R_0 = R0, pop_0 = pop0) %>%
      select(parid, rho, tau, beta1, beta2, beta3, beta4, beta5,
                    beta6, nu, gamma, sigma, theta0, alpha, mu, delta, sig_sq,
                    S_0, E_0, I_0, A_0, R_0, pop_0)

    write_csv(starts_epi, paste0(out_dir, "starts_", mcode, "_epi", nweeks,
                                 paste0("_nprofs", nprofs), "_", parseed, ".csv"))
  }
  return(starts_epi)
}

#### set parameters for endemic fit model
generate_end_params <- function(out_dir_epi, out_dir_end, episettings, endsettings) {
  set.seed(endsettings$seed)
  est_epi <- readRDS(paste0(out_dir_epi, "iffits_epi.rds"))
  if(file.exists(paste0(out_dir_end, "starts_", endsettings$mcode, "_epi",
                        episettings$nweeks, ".csv"))) {
    starts_end <- read_csv(paste0(out_dir_end, "starts_", endsettings$mcode,
                                  "_epi", episettings$nweeks, ".csv"))

  }
  else {
    starts_end1 <- est_epi %>%
      dplyr::mutate(epi_mcode = episettings$mcode) %>%
      dplyr::mutate(parid = as.integer(parid)) %>%
      dplyr::select(epi_mcode, parid, loglik, rho, tau, beta1, beta2, beta3,
                    beta4, beta5, beta6, nu, gamma, sigma, theta0, alpha, mu, sig_sq,
                    delta)
    fitstate_fns <- list.files(out_dir_epi, paste0("fitstates_if_",
                                                   episettings$str_if2))

    initconds <- map_dfr(1:length(fitstate_fns), function(i) {
      par_id <- gsub(paste0("_seed", episettings$seed, ".csv"), "",
                     gsub(paste0("fitstates_if_",
                                 episettings$str_if2, "_parid"),
                          "", fitstate_fns[i]))
      read_csv(paste0(out_dir_epi, fitstate_fns[i])) %>%
        dplyr::filter(week == max(week)) %>%
        dplyr::mutate(S_0 = S_med/pop_med, E_0 = E_med/pop_med,
                      I_0 = I_med/pop_med, A_0 = A_med/pop_med,
                      incid_0 = incid_med) %>%
        dplyr::mutate(R_0 = 1-(S_0+E_0+I_0+A_0), parid = as.numeric(par_id),
                      pop_0 = pop_med) %>%
        dplyr::select(parid, S_0, E_0, I_0, A_0, R_0, incid_0, pop_0)
    }) %>% arrange(parid)

    starts_end <- inner_join(starts_end1, initconds, by = c("parid")) %>%
      dplyr::filter(!is.na(S_0)) %>%
      dplyr::filter(nu > 0.9 & beta1 < 100) %>%
      dplyr::arrange(desc(loglik))
    write_csv(starts_end, paste0(out_dir_end, "starts_", endsettings$mcode,
                                 "_epi", episettings$nweeks, ".csv"))
  }
  return(starts_end)
}

##############################################################################--
####  UTILITY FUNCTIONS FOR MODEL ASSESSMENT                                ####
##############################################################################--
#### simulate data and summarise over simulations
sim_data <- function(sub_model, pars, num_sims) {
  sims <- simulate(sub_model, params = pars,
                   nsim = num_sims, format = "data.frame") %>%
    dplyr::mutate(pop = S + E + I + A + R) %>%
    dplyr::select(week, S, E, I, A, R, incid, pop, cases)

  states <- sims %>%
    dplyr::group_by(week) %>%
    dplyr::summarise(S_med = median(S), E_med = median(E), I_med = median(I),
                     A_med = median(A), R_med = median(R),
                     incid_med = median(incid), pop_med = median(pop),
                     cases_med = median(cases),
                     S_mean = mean(S), E_mean = mean(E), I_mean = mean(I),
                     A_mean = mean(A), R_mean = mean(R),
                     incid_mean = mean(incid), pop_mean = mean(pop),
                     cases_mean = mean(cases),
                     S_lo = quantile(S, probs = c(.025)),
                     E_lo = quantile(E, probs = c(.025)),
                     I_lo = quantile(I, probs = c(.025)),
                     A_lo = quantile(A, probs = c(.025)),
                     R_lo = quantile(R, probs = c(.025)),
                     incid_lo = quantile(incid, probs = c(.025)),
                     pop_lo = quantile(pop, probs = c(.025)),
                     cases_lo = quantile(cases, probs = c(.025)),
                     S_hi = quantile(S, probs = c(.975)),
                     E_hi = quantile(E, probs = c(.975)),
                     I_hi = quantile(I, probs = c(.975)),
                     A_hi = quantile(A, probs = c(.975)),
                     R_hi = quantile(R, probs = c(.975)),
                     incid_hi = quantile(incid, probs = c(.975)),
                     pop_hi = quantile(pop, probs = c(.975)),
                     cases_hi = quantile(cases, probs = c(.975))) %>%
    dplyr::mutate(date = lubridate::ymd("2010-10-14") + lubridate::weeks(week))
  return(states)
}

#### record results as min, max, med, mean for estimated parameters and states
rec_res <- function(res_df) {
  val_ests_all <- res_df %>%
    dplyr::select(parid, rho, tau, beta1, beta2, beta3, beta4, beta5, beta6,
                  nu, S_0, E_0, I_0, A_0, R_0, loglik) %>%
    dplyr::mutate(S_0 = round(S_0 * pop.haiti),
                  E_0 = round(E_0 * pop.haiti),
                  I_0 = round(I_0 * pop.haiti),
                  A_0 = round(A_0 * pop.haiti),
                  R_0 = round(R_0 * pop.haiti))
  medians <- apply(val_ests_all, 2, median) %>% round(2)
  means <- apply(val_ests_all, 2, mean) %>% round(2)
  maxima <- apply(val_ests_all, 2, max) %>% round(2)
  minima <- apply(val_ests_all, 2, min) %>% round(2)
  val_ests <- data.frame(colnames(val_ests_all), minima, maxima, medians, means)
  colnames(val_ests) <- c("parameter", "minimum", "maximum", "median", "mean")
  return(val_ests)
}

#### read in simulated data for fitted states and calculate quantiles
process_fitstates <- function(fitstate_fns, hti_data) {
  map_dfr(1:length(fitstate_fns), function(i){
    read_csv(fitstate_fns[i], col_types = cols())
  }) %>%
    ungroup %>%
    group_by(week) %>%
    dplyr::summarise(est_med = median(cases_med),
                     est_lo = quantile(cases_med, probs = c(.025)),
                     est_hi = quantile(cases_med, probs = c(.975))) %>%
    dplyr::full_join(hti_data %>% dplyr::rename(orig_cases = cases),
                     by = c("week") )%>%
    dplyr::arrange(week) -> t
  t$date <- lubridate::ymd("2010-10-23") + lubridate::weeks(t$week)
  return(t)
}
