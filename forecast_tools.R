#==============================================================================#
# Model 1: Helper functions for forecasting using model 1                      #
#          Code modified from Elizabeth C. Lee, Andrew S. Azman, and Justin    #
#          Lessler at Johns Hopkins Bloomberg School of Public Health          #
#==============================================================================#

##----model projection settings-----------------------------------------------##
fix_fc_settings <- function(nsims) {
  modcode <- "final"
  scencode <- "novac"
  fc_seed <- 20192104
  fc_nsims <- nsims
  fc_horizon <- 52*11 ## 11 years forecast
  
  
  fc_str_if <- paste0(scencode, "_nsims", fc_nsims, "_seed", fc_seed)
  fc_str_if2 <- paste0(scencode, "_nsims", fc_nsims)
  
  return(list(mcode=modcode, scode=scencode, seed=fc_seed, 
              horizon = fc_horizon, nsims=fc_nsims, str_if = fc_str_if, 
              str_if2 = fc_str_if2))
}

fix_fc_vac_settings <- function(scen_code, nsims) {
  modcode <- "final"
  scencode <- scen_code
  fc_seed <- 20192104
  fc_nsims <- nsims
  fc_horizon <- 52*11 ## 11 years forecast
  
  fc_nd = NA; fc_nw = NA; fc_c2 = NA; fc_c1 = NA; fc_vescen = NA
  
  ## grouped by deployment strategy
  if (scencode %in% paste0("id", seq(2, 34, by = 4))){ ## 2 dept
    fc_nd = 2; fc_nw = 52
  } else if (scencode %in% paste0("id", seq(4, 36, by = 4))){ ## 3 dept 
    fc_nd = 3; fc_nw = 33
  } else if (scencode %in% paste0("id", seq(3, 35, by = 4))){ ## slow national
    fc_nd = 10; fc_nw = 26
  } else if (scencode %in% paste0("id", seq(1, 33, by = 4))){ ## fast national 
    fc_nd = 10; fc_nw = 10
  } else {
    fc_nd = 0; fc_nw = 26
  }
  
  ## grouped by coverage level
  if (scencode %in% paste0("id", 1:12)){
    fc_c2 = 0.7; fc_c1 = 0.1
  } else if (scencode %in% paste0("id", 13:24)){
    fc_c2 = 0.4; fc_c1 = 0.2
  } else if (scencode %in% paste0("id", 25:36)){
    fc_c2 = 0.95; fc_c1 = 0.0167
  } else {
    fc_c2 = 0; fc_c1 = 0
  }
  
  ## grouped by ve scenario
  if (scencode %in% paste0("id", c(1:4, 13:16, 25:28))){
    fc_vescen = "ve_s1"
  } else if (scencode %in% paste0("id", c(5:8, 17:20, 29:32))){
    fc_vescen = "ve_s2"
  } else if (scencode %in% paste0("id", c(9:12, 21:24, 33:36))){
    fc_vescen = "ve_s3"
  }
  
  fc_str_if <- paste0(scencode, "_nsims", fc_nsims, "_seed", fc_seed)
  fc_str_if2 <- paste0(scencode, "_nsims", fc_nsims)
  
  return(list(mcode = modcode, scode = scencode, seed = fc_seed, horizon = fc_horizon, 
              nsims = fc_nsims, nd = fc_nd, nw = fc_nw, c2 = fc_c2, c1 = fc_c1, 
              vescen = fc_vescen, str_if = fc_str_if, str_if2 = fc_str_if2))
}

##----set parameters for projection model-------------------------------------##
generate_fc_params <- function(out_dir_epi, out_dir_end, out_dir_fc,
                               episettings, endsettings, fcsettings) {
  set.seed(fcsettings$seed)
  
  kap <- 0.95
  depts <- fcsettings$nd
  
  est_epi <- readRDS(paste0(out_dir_epi, "iffits_", 
                            episettings$str_if, ".rds")) %>%
    dplyr::select(parid, loglik) %>%
    dplyr::filter(!is.na(loglik)) %>%
    dplyr::filter(is.finite(loglik)) %>%
    dplyr::rename(loglik_epi = loglik)
  est_end <- readRDS(paste0(out_dir_end, "iffits_", 
                            endsettings$str_if, ".rds")) %>%
    dplyr::filter(!is.na(loglik)) %>%
    dplyr::filter(is.finite(loglik)) %>%
    dplyr::mutate(mcode = endsettings$mcode) %>%
    dplyr::mutate(parid = as.integer(parid)) %>%
    left_join(est_epi, by = c("parid")) %>%
    dplyr::mutate(loglik_f = loglik + loglik_epi) %>%
    dplyr::select(mcode, parid, loglik_f, rho, tau, beta1, beta2, beta3, 
                  beta4, beta5, beta6, nu, gamma, sigma, theta0, alpha, 
                  mu, delta) %>%
    dplyr::mutate(kappa = kap)
  
  incl_parids <- unlist(est_end$parid)
  fitstate_fns <- paste0("fitstates_if_", endsettings$str_if2, "_parid", 
                         incl_parids, "_seed", endsettings$seed, ".csv")
  fitstate_fns <- unique(fitstate_fns)
  initconds <- map_dfr(1:length(fitstate_fns), function(i){
    par_id <- incl_parids[i]
    read_csv(paste0(out_dir_end, fitstate_fns[i]), col_types = cols()) %>%
      dplyr::filter(week == max(week)) %>%
      dplyr::mutate(S_0 = S_med/pop_med, E_0 = E_med/pop_med, 
                    I_0 = I_med/pop_med, A_0 = A_med/pop_med, 
                    incid_0 = incid_med) %>%
      dplyr::mutate(R_0 = 1-(S_0+E_0+I_0+A_0), parid = as.numeric(par_id), 
                    pop_0 = pop_med) %>%
      dplyr::mutate(E_0 = ifelse(E_0 == 0, 1E-7, E_0), 
                    I_0 = ifelse(I_0 == 0, 1E-7, I_0)) %>%
      dplyr::select(parid, S_0, E_0, I_0, A_0, R_0, incid_0, pop_0)
  })
  
  starts_fc <- inner_join(est_end, initconds, by = c("parid")) %>%
    dplyr::filter(!is.na(S_0)) %>%
    dplyr::arrange(desc(loglik_f))
  write_csv(starts_fc, paste0(out_dir_fc, "starts_", fcsettings$mcode, 
                              "_epi", episettings$nweeks, ".csv"))
  if (depts > 1) {
    init_states <- c(paste0("S", 1:depts, "_0"), paste0("E", 1:depts, "_0"),
                     paste0("I", 1:depts, "_0"), paste0("A", 1:depts, "_0"),
                     paste0("R", 1:depts, "_0"))
    for (i in 1:length(init_states)) {
      columns <- colnames(starts_fc)
      columns <- c(columns, init_states[i])
      starts_fc[, (ncol(starts_fc) + 1)] <- 0.0
      colnames(starts_fc) <- columns
    }
  }
  
  return(starts_fc)
}


#### construct vaccination covariate table
make_vactab <- function(t0 = 0, tmax, byt = 1, ndept = 10, nweeks = 10, 
                        coverage_2dose, coverage_1dose, first_vac_t, ve_scen) {
  
  ## deployment scenario 1 (fast national): one departmental campaign every 10 weeks for all depts
  
  tbasis <- seq(from=t0,to=tmax,by=byt)
  time_check <- c()
  for(i in 1:ndept) {
    time_check <- c(time_check, rep(0, nweeks-1), rep(i, 1))
  } ## 10 week vacc campaigns (nweeks)
  
  ## number of vaccines per week by department
  pop_dept <- data.frame(ocv_order = 1:10, 
                         dept = c("Centre", "Artibonite", "Ouest", "Nord Ouest", 
                                  "Nord", "Sud", "Nippes", "Nord Est", "Sud Est", 
                                  "Grand'Anse"), 
                         pop = c(746236, 1727524, 4029705, 728807, 1067177, 
                                 774976, 342525, 393967, 632601, 468301)) %>%
    mutate(num_vacc = (coverage_2dose+coverage_1dose)*pop/1) ## pulse vaccinees in last 1 week of campaign
  # browser()
  ## create dataframe with number vaccinated for each campaign
  vactab <- data.frame(time = tbasis, vac_tcheck = 0)
  vactab[which(vactab$time %in% first_vac_t:(first_vac_t+length(time_check)-1)),]$vac_tcheck <- time_check
  vactab2 <- left_join(vactab, pop_dept %>% select(ocv_order, num_vacc), by = c("vac_tcheck"="ocv_order")) %>%
    mutate(num_vacc = ifelse(is.na(num_vacc), 0, round(num_vacc))) 
  
  ## ve decay after x weeks
  veDecay_mo <- read_csv("Data/ve_decay_bymonth.csv", col_types = cols())
  veDecay_adult <- bind_rows(veDecay_mo, veDecay_mo, veDecay_mo, veDecay_mo) %>% 
    arrange(month) %>%
    select(!!ve_scen) %>% 
    unlist %>% unname
  ## adjust for lower VE in U5 population (U5 VE is 0.4688*adult VE; roughly 11% of the population is 0-4 years old according to UN World Population Prospects 2017)
  ## adjust ve for age (population VE = adult VE * (1-(1-0.4688)*proportion under-5))
  veDecay <- veDecay_adult * (1-(1-0.4688)*0.11)
  ## adjust for one-dose decay after 52 weeks
  veDecay[53:length(veDecay)] <- coverage_2dose/(coverage_2dose+coverage_1dose)*veDecay[53:length(veDecay)]
  
  ## add vaccine immunity decay, time shifted for each campaign
  decay_times <- vactab2 %>%
    filter(vac_tcheck > 0) %>%
    group_by(vac_tcheck) %>%
    filter(time == max(time)) %>% 
    select(-num_vacc) %>%
    ungroup
  for(i in 1:ndept) {
    decay_start <- decay_times %>%
      filter(vac_tcheck == i) %>% 
      select(time) %>% unlist %>% unname
    decay_df <- as_tibble(data.frame(time = (seq_along(veDecay)+decay_start))) %>%
      mutate(!!paste0("ve_d", i) := veDecay)
    vactab2 <- left_join(vactab2, decay_df, by = c("time")) 
  }
  
  vactab3 <- vactab2 %>%
    mutate_at(vars(contains("ve_")), list(~ifelse(is.na(.), 0, .)))
  
  return(vactab3)
}

## -------------------------------------- ##
## process data for elimination
process_fc_elim <- function(fcsettings, fcstate_fns) {
  ## fcstate_fns include path
  ## sim, week, cases, S, E, I, A, R, N, incid, loglik
  map_dfr(1:length(fcstate_fns), function(i){
    par_id <- gsub(paste0("_seed", fcsettings$seed, ".csv"), "", 
                   gsub(paste0("=/GeneratedData/", fcsettings$mcode, "_", 
                               fcsettings$scode, "/fcstates_stoch_", 
                               fcsettings$str_if2, "_parid"), "", fcstate_fns[i]))
    print(paste("***** par id", par_id, "*****"))
    if (is.na(fcsettings$vescen)) { ## no vaccinations
      dummy <- read_csv(fcstate_fns[i], col_types = cols()) %>%
        dplyr::select(sim, week, cases, incid, pop)
    } else { ## yes vaccinations
      dummy <- read_csv(fcstate_fns[i], col_types = cols()) %>%
        dplyr::select(sim, week, cases, incid, asymV, pop)
    }
    output <- find_elimPeriod(fcsettings, dummy) %>%
      dplyr::mutate(parid = as.numeric(par_id))
    rm(dummy)
    gc()
    return(output)
  }) 
}

find_elimPeriod <- function(fcsettings, simdata) {
  vacstart <- min(simdata$week) + 4
  ## should work with sims_tm and sims_if
  simdata2 <- simdata %>% 
    dplyr::filter(week >= vacstart)
  if (is.na(fcsettings$vescen)) { ## no vaccinations
    simdata2 <- simdata2 %>%
      dplyr::mutate(true_incid = incid) %>%
      ungroup 
  } else { ## yes vaccinations
    simdata2 <- simdata2 %>%
      dplyr::mutate(true_incid = incid + asymV) %>%
      ungroup 
  }
  uqsims <- unique(simdata2$sim)
  
  returnData <- map_dfr(uqsims, function(uqsim) {
    rleSimdata <- simdata2 %>% dplyr::filter(sim == uqsim) %>%
      dplyr::mutate(window_marker = seq_along(sim))
    print(paste("uqsim", uqsim, "**********"))
    
    ## 1/10000 incidence threshold
    rleinit <- rleSimdata %>% dplyr::filter(window_marker %in% 1:52)
    while (sum(rleinit$true_incid)>=(1/10000*mean(rleinit$pop)) & nrow(rleinit)>=52) {
      ## use the while loop to find and filter the minimum (52 weeks) elimination dataset
      new_window_start <- min(rleinit$window_marker)+1
      new_window_end <- new_window_start+51
      rleinit <- rleSimdata %>% 
        dplyr::filter(window_marker %in% new_window_start:new_window_end)
    }
    if (nrow(rleinit) < 52) { ## no elimination period dataset
      dummyA <- rep(F, nrow(rleSimdata))
      resA <- NA
      print(paste("no elim @ 1/10K threshold"))
    } else { ## find true elimination period dataset
      elim_testresurg_start <- max(rleinit$window_marker)+1
      testresurg <- rleSimdata %>% 
        dplyr::filter(window_marker %in% elim_testresurg_start:(elim_testresurg_start+52))
      dummyA <- rep(F, nrow(rleSimdata))
      resA <- FALSE
      print(paste("elim achieved @ 1/10K threshold"))
      while (sum(testresurg$true_incid)<(1/10000*mean(testresurg$pop)) & max(testresurg$window_marker)<max(rleSimdata$window_marker)) {
        new_window_start <- min(testresurg$window_marker)+1
        new_window_end <- new_window_start+51
        testresurg <- rleSimdata %>% 
          dplyr::filter(window_marker %in% new_window_start:new_window_end)
      }
      if (sum(testresurg$true_incid)>=(1/10000*mean(testresurg$pop))) { ## resurgence occurs
        ## find entire duration of elimination period (prior to resurgence)
        elim_beforeresurg_start <- min(rleinit$window_marker)
        elim_beforeresurg_end <- min(testresurg$window_marker)-1
        beforeresurg <- rleSimdata %>% 
          dplyr::filter(window_marker %in% elim_beforeresurg_start:elim_beforeresurg_end)
        dummyA[which(rleSimdata$window_marker %in% beforeresurg$window_marker)] <- TRUE
        resA <- TRUE
        print(paste("resurgence occurs"))
      } else { ## resurgence does not occur
        elim_start <- min(rleinit$window_marker)
        elim_end <- max(testresurg$window_marker)
        noresurg <- rleSimdata %>% 
          dplyr::filter(window_marker %in% elim_start:elim_end) 
        dummyA[which(rleSimdata$window_marker %in% noresurg$window_marker)] <- TRUE
        resA <- FALSE
        print(paste("resurgence does not occur"))
      }
    }
    
    ## 1 case threshold
    print("processing one case threshold")
    dummyB <- rep(F, nrow(rleSimdata))
    resB <- NA
    rle.resultsB <- rle(rleSimdata$true_incid)
    if (rle.resultsB$values[length(rle.resultsB$values)]==0 & rle.resultsB$lengths[length(rle.resultsB$values)] > 52) { 
      ## elimination period that is at least 52 weeks long and endures through the end of the 10 year simulation
      first.rle.indexB = length(rle.resultsB$values)
      pre.indexB = ifelse(first.rle.indexB>1, sum(rle.resultsB$lengths[1:(first.rle.indexB-1)])+1, 1)
      post.indexB = pre.indexB + rle.resultsB$lengths[first.rle.indexB]-1
      dummyB[pre.indexB:post.indexB] <- T
      resB <- NA ## 6/11/2019 resurgence is not possible in this new definition of elimination
      print(paste("elim @ 1 case threshold in sim", uqsim))
    } else {
      print(paste("no elim @ 1 case threshold"))
    }
    
    rleSimdata %>% 
      dplyr::mutate(elimPeriodA = dummyA, elimPeriodB = dummyB, resurgA = resA, 
                    resurgB = resB) %>% ## is there resurgence to non-elimination levels for a given simulation (same value for all time steps)
      dplyr::select(sim, week, true_incid, pop, elimPeriodA, elimPeriodB, resurgA, resurgB)
  })

  return(returnData)
}

find_probElim <- function(elimdata) {
  year3 <- round(min(elimdata$week) + 3*52.14)
  year5 <- round(min(elimdata$week) + 5*52.14)
  year10 <- round(min(elimdata$week) + 10*52.14)
  
  y3 <- elimdata %>%
    dplyr::filter(week == year3) %>%
    summarise(nElimA = sum(elimPeriodA), 
              nElimB = sum(elimPeriodB), 
              nsims = nrow(.), 
              resurgA = sum(resurgA, na.rm=TRUE),
              resurgB = sum(resurgB, na.rm=TRUE)) %>%
    dplyr::mutate(pr_elim_thresh1 = nElimA/nsims, 
                  pr_elim_thresh2 = nElimB/nsims,
                  pr_resurg_thresh1 = resurgA/nElimA,
                  pr_resurg_thresh2 = resurgB/nElimB) %>%
    dplyr::mutate(time_frame = 3)
  y5 <- elimdata %>%
    dplyr::filter(week == year5) %>%
    summarise(nElimA = sum(elimPeriodA), 
              nElimB = sum(elimPeriodB), 
              nsims = nrow(.), 
              resurgA = sum(resurgA, na.rm=TRUE),
              resurgB = sum(resurgB, na.rm=TRUE)) %>%
    dplyr::mutate(pr_elim_thresh1 = nElimA/nsims, 
                  pr_elim_thresh2 = nElimB/nsims, 
                  pr_resurg_thresh1 = resurgA/nElimA, 
                  pr_resurg_thresh2 = resurgB/nElimB) %>%
    dplyr::mutate(time_frame = 5) 
  y10 <- elimdata %>%
    dplyr::filter(week == year10) %>%
    summarise(nElimA = sum(elimPeriodA), 
              nElimB = sum(elimPeriodB), 
              nsims = nrow(.), 
              resurgA = sum(resurgA, na.rm=TRUE),
              resurgB = sum(resurgB, na.rm=TRUE)) %>%
    dplyr::mutate(pr_elim_thresh1 = nElimA/nsims, 
                  pr_elim_thresh2 = nElimB/nsims, 
                  pr_resurg_thresh1 = resurgA/nElimA, 
                  pr_resurg_thresh2 = resurgB/nElimB) %>%
    dplyr::mutate(time_frame = 10) 
  
  bind_rows(y3, y5, y10)
}

find_timeToElim <- function(elimdata) {
  vacstart <- min(elimdata$week)
  periodA <- elimdata %>%
    dplyr::filter(elimPeriodA) %>%
    group_by(parid, sim) %>% 
    dplyr::filter(week == min(week)) %>%
    ungroup %>%
    dplyr::mutate(tElim_wks = week-vacstart) %>%
    summarise(t_elim_thresh1_med = median(tElim_wks),
              t_elim_thresh1_low = quantile(tElim_wks, probs = c(.025)), 
              t_elim_thresh1_hi = quantile(tElim_wks, probs = c(.975)))
  
  periodB <- elimdata %>%
    dplyr::filter(elimPeriodB) %>%
    group_by(parid, sim) %>% 
    dplyr::filter(week == min(week)) %>%
    ungroup %>%
    dplyr::mutate(tElim_wks = week-vacstart) %>%
    summarise(t_elim_thresh2_med = median(tElim_wks), 
              t_elim_thresh2_low = quantile(tElim_wks, probs = c(.025)), 
              t_elim_thresh2_hi = quantile(tElim_wks, probs = c(.975))) 
  
  bind_cols(periodA, periodB)
}

################################################




## -------------------------------------- ##
## process data for plotting

process_fc_plot <- function(fitstate_fns, hti_data) {
  ## fitstate_fns include path
  map_dfr(1:length(fitstate_fns), function(i){
    read_csv(fitstate_fns[i], col_types = cols()) 
  }) %>%
    ungroup %>%
    group_by(week) %>%
    dplyr::rename(cmed = cases_med, imed = incid_med) %>%
    dplyr::summarise(cases_med = median(cmed), 
                     cases_lo = quantile(cmed, probs = c(.025)), 
                     cases_hi = quantile(cmed, probs = c(.975)), 
                     incid_med = median(imed), 
                     incid_lo = quantile(imed, probs = c(.025)), 
                     incid_hi = quantile(imed, probs = c(.975)), 
                     pop_med = median(pop_med)) %>%
    dplyr::full_join(hti_data %>% dplyr::rename(orig_cases = cases), by = c("week"))
}

process_ffc_plot <- function(fcstate_fns){
  ## fcstate_fns include path
  ## sim, week, cases, S, E, I, A, R, N, incid, loglik
  map_dfr(1:length(fcstate_fns), function(i){
    readRDS(fcstate_fns[i])
    # read_csv(fcstate_fns[i])
  }) %>%
    dplyr::mutate(orig_cases = NA) %>%
    dplyr::select(week, cases_med, cases_lo, cases_hi, incid_med, incid_lo, incid_hi, orig_cases, pop_med)
}


## -------------------------------------- ##
## process cumulative vaccination plot

get_vac_settings <- function(scencode, ndata){
  ## grouped by deployment strategy
  if (scencode %in% paste0("id", seq(2, 34, by = 4))){ ## 2 dept campaigns
    fc_nd = 2; fc_nw = 52
  } else if (scencode %in% paste0("id", seq(4, 36, by = 4))){ ## 3 dept campaigns
    fc_nd = 3; fc_nw = 33
  } else if (scencode %in% paste0("id", seq(3, 35, by = 4))){ ## slow national campaigns
    fc_nd = 10; fc_nw = 26
  } else if (scencode %in% paste0("id", seq(1, 33, by = 4))){ ## fast national campaigns
    fc_nd = 10; fc_nw = 10
  } 
  
  ## grouped by coverage level
  if (scencode %in% paste0("id", 1:12)){
    fc_c2 = 0.7; fc_c1 = 0.1
  } else if (scencode %in% paste0("id", 13:24)){
    fc_c2 = 0.4; fc_c1 = 0.2
  } else if (scencode %in% paste0("id", 25:36)){
    fc_c2 = 0.95; fc_c1 = 0.0167
  }
  
  ## grouped by ve scenario
  if (scencode %in% paste0("id", c(1:4, 13:16, 25:28))){
    fc_vescen = "ve_s1"
  } else if (scencode %in% paste0("id", c(5:8, 17:20, 29:32))){
    fc_vescen = "ve_s2"
  } else if (scencode %in% paste0("id", c(9:12, 21:24, 33:36))){
    fc_vescen = "ve_s3"
  }
  
  vactab <- make.vactab(t0 = 0, 
                        tmax = ndata+53*6, 
                        ndept = fc_nd,
                        nweeks = fc_nw,
                        coverage_2dose = fc_c2, 
                        coverage_1dose = fc_c1,
                        first_vac_t = ndata+4, 
                        ve_scen = fc_vescen) %>%
    dplyr::select(time, num_vacc) %>% 
    dplyr::filter(time > ndata)
  
  return(vactab)
}

get_date_df <- function(){
  ndp <- nrow(get.mspp.agg.data())
  
  date.start <- get.mspp.dept.data() %>% dplyr::filter(date_sat == min(date_sat)) %>% dplyr::distinct(date_sat)
  datedf <- data.frame(saturday_date = seq(date.start$date_sat, length.out = ndp+52*11, by = 7), week = seq_len(ndp+52*11))
  return(datedf)
}

process_cumvac_output <- function(){
  ndp <- nrow(get.mspp.agg.data())
  
  date.start <- get.mspp.dept.data() %>% dplyr::filter(date_sat == min(date_sat)) %>% dplyr::distinct(date_sat)
  datedf <- data.frame(weekdate = seq(date.start$date_sat, length.out = ndp+52*11, by = 7), week = seq_len(ndp+52*11))
  
  s1v <- get_vac_settings("id1", ndp) %>% dplyr::mutate(scenario = 1)
  s2v <- get_vac_settings("id2", ndp) %>% dplyr::mutate(scenario = 2)
  s3v <- get_vac_settings("id3", ndp) %>% dplyr::mutate(scenario = 3)
  s4v <- get_vac_settings("id4", ndp) %>% dplyr::mutate(scenario = 4)
  
  vacc_scens <- bind_rows(s1v, s2v, s3v, s4v) %>%
    rename(week = time) %>%
    dplyr::mutate(scenario = as.character(scenario)) %>%
    left_join(datedf, by = c("week")) %>%
    dplyr::rename(vacc_med = num_vacc)
  
  return(vacc_scens)
}

process_cumvac_plot <- function(datedf){
  ndp <- nrow(get.mspp.agg.data())
  s1v <- get_vac_settings("id1", ndp) %>% dplyr::rename(s1 = num_vacc)
  s2v <- get_vac_settings("id2", ndp) %>% dplyr::rename(s2 = num_vacc)
  s3v <- get_vac_settings("id3", ndp) %>% dplyr::rename(s3 = num_vacc)
  s4v <- get_vac_settings("id4", ndp) %>% dplyr::rename(s4 = num_vacc)
  
  vacc_scens <- full_join(s1v, s2v, by = c("time")) %>%
    full_join(s3v, by = c("time")) %>%
    full_join(s4v, by = c("time")) %>%
    rename(week = time) %>%
    left_join(datedf, by = c("week")) 
  
  return(vacc_scens)
}