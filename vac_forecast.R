#==============================================================================#
# Model 1: Vaccination Forecasting                                             #
#          Predicting resurgence and elimination with a vaccination campaign   #      
#          (2 dept., 3 dept., fast nat., slow nat., slow 2 dept.)              #
#          Code modified from Elizabeth C. Lee, Andrew S. Azman, and Justin    #
#          Lessler at Johns Hopkins Bloomberg School of Public Health          #
#==============================================================================#
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
find_rtools()
has_devel()
source("mod_build.R")
source("mod_fc_build.R")
source("forecast_tools.R")
source("mod_fit.R")
source("mod_vis.R")
source("mod_tools.R")
registerDoParallel(cores=6)

##----settings --- change if desired--------------------------------------------
# filtering and simulation settings
num_iters <- 50 # number of iterations for IF2 (originally 100)
num_parts <- 1000 # number of particles for IF2 (originally 100)
num_sims <- 25 # number of simulations for median calculations (originally 25)

# simulation evaluation settings
plot_modules <- TRUE # plot 
process_modules <- TRUE
cumvac_modules <- TRUE

# overall evaluation settings
ts_module <- TRUE
pr_elim_module <- TRUE # find 
t_elim_module <- TRUE # find time to elimination 
numvacc_module <- TRUE # find approximate number of doses needed

# make directory for data generated from simulations and filtering
dir.create("GeneratedData", showWarnings = FALSE)

# make directory for figures
dir.create("figures", showWarnings = FALSE)

# make directory for csv output
dir.create("GeneratedData/final1_summaries", showWarnings = FALSE)

## NOTE: # represents the code of a different vaccination deployment strategy 
## (1 = fast national, 2 = 2-department, 3 = slow national, 4 = 3-department, 
##  25 = fast national high coverage, and all others represent combinations)
scen_code = "id2"

####VAC MODEL FORECASTING#######################################################
episettings <- fix_epi_settings(nprofs = 3, nparticles = 1000, nmif = 150)
endsettings <- fix_end_settings(episettings)
fcsettings <- fix_fc_vac_settings(scen_code = "id2", nsims = 25)
scode <- fcsettings$scode

fig_dir <- paste0("figures/", fcsettings$mcode, "_", scode, "/")
out_dir_fc <- paste0("GeneratedData/", fcsettings$mcode, "_", scode, "/")
out_dir_end <- paste0("GeneratedData/", endsettings$mcode, "/")
out_dir_epi <- paste0("GeneratedData/", episettings$mcode, "/")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir_fc, showWarnings = FALSE)

##----import settings and data------------------------------------------------##
starts <- generate_fc_params(out_dir_epi, out_dir_end, out_dir_fc,
                             episettings, endsettings, fcsettings)
pop.haiti <- get_haiti_pop()
haiti.dat <- get_mspp_agg_data() 
haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= episettings$nweeks)
haiti.dat.end <- haiti.dat %>% dplyr::filter(week > episettings$nweeks)

depts <- fcsettings$nd

allpars <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6",
             "gamma","sigma","theta0","alpha","mu","delta","nu","kappa",
             "S_0","E_0","I_0","A_0","R_0",
             paste0("S", 1:depts, "_0"),
             paste0("E", 1:depts, "_0"),
             paste0("I", 1:depts, "_0"),
             paste0("A", 1:depts, "_0"),
             paste0("R", 1:depts, "_0"),
             "incid_0","pop_0")

## vac covariates
vactab <- make_vactab(t0 = 0, 
                      tmax = nrow(haiti.dat) + fcsettings$horizon + 1, 
                      ndept = fcsettings$nd,
                      nweeks = fcsettings$nw,
                      coverage_2dose = fcsettings$c2, 
                      coverage_1dose = fcsettings$c1,
                      first_vac_t = nrow(haiti.dat) + 4, 
                      ve_scen = fcsettings$vescen)

## spline covariates
covartab.all <- make_covartab(0, nrow(haiti.dat) + fcsettings$horizon + 1, 
                              byt = 1, degree = 6, nbasis = 6, per = 52.14) %>%
  dplyr::mutate(time_endfc = ifelse(time < episettings$nweeks, NA, time)) %>%
  dplyr::mutate(time_fc = ifelse(time < max(haiti.dat$week), NA, time)) %>%
  dplyr::mutate(time_end = ifelse(time <= max(haiti.dat$week), time, NA))
covartab.all2 <- covartab.all %>%
  dplyr::select(time, contains("seas")) %>%
  full_join(vactab, by = c("time"))

##----make model----------------------------------------------------------------
start_pars <- starts %>% dplyr::select(-mcode)
haiti.mod.fc <- mod_fc_build(dat = haiti.dat,
                             my.times = "week",
                             covar.times = "time",
                             my.t0 = 0, 
                             covar = covartab.all2, 
                             depts = depts)

##----perform stochastic forecasts----------------------------------------------
bake(file=paste0(out_dir_fc, "/stochfcs_", fcsettings$str_if, ".rds"), 
     seed = fcsettings$seed, {
  foreach(start = iter(start_pars, by = 'row'),
          .inorder = FALSE, .combine = rbind,
          .packages = c("pomp", "magrittr"),
          .errorhandling = c('remove'),
          .noexport = c(),
          .verbose = TRUE) %dopar%
    {
      tic <- Sys.time()
      print(paste("forecast stoch params vac", start$parid))
      
      M2 <- haiti.mod.fc
      time(M2) <- seq((tail(haiti.dat.end$week, n = 1) + 1), tail(covartab.all2$time, n = 1))
      timezero(M2) <- tail(haiti.dat.end$week, n = 1)
      sim_pars <- unlist(start[which(names(start) %in% allpars)])
    
      base_sims <- simulate(M2, nsim = fcsettings$nsims, params = sim_pars, format = "data.frame") 
      base_sims <- simulate(M2, params = sim_pars, format = "data.frame") 
      
      if (scode %in% paste0("id", seq(2, 34, by = 4))) { ## 2 dept
        fc_states_sims <- base_sims %>%
          dplyr::mutate(pop_nv = S+E+I+A+R, 
                        pop_1 = S1+E1+I1+A1+R1, 
                        pop_2 = S2+E2+I2+A2+R2, 
                        pop = pop_nv+pop_1+pop_2, 
                        loglik = start$loglik_f) %>%
          #dplyr::rename(week = time) %>%
          dplyr::rename(sim = .id) %>%
          dplyr::select(sim, week, cases, 
                        S, E, I, A, R, pop_nv, 
                        S1, E1, I1, A1, R1, pop_1, 
                        S2, E2, I2, A2, R2, pop_2, 
                        pop, incid, incidU, incidV, asymV, newV, loglik, 
                        foival, Str0, Sout, Sin)
      } else if (scode %in% paste0("id", seq(4, 36, by = 4))) { ## 3 dept
        fc_states_sims <- base_sims %>%
          dplyr::mutate(pop_nv = S+E+I+A+R, 
                        pop_1 = S1+E1+I1+A1+R1, 
                        pop_2 = S2+E2+I2+A2+R2, 
                        pop_3 = S3+E3+I3+A3+R3,
                        pop = pop_nv+pop_1+pop_2+pop_3, 
                        loglik = start$loglik_f) %>%
          # dplyr::rename(week = time) %>%
          dplyr::rename(sim = .id) %>%
          dplyr::select(sim, week, cases, 
                        S, E, I, A, R, pop_nv, 
                        S1, E1, I1, A1, R1, pop_1, 
                        S2, E2, I2, A2, R2, pop_2,
                        S3, E3, I3, A3, R3, pop_3, 
                        pop, incid, incidU, incidV, asymV, newV, loglik, 
                        foival, Str0, Sout, Sin)
      } else if (scode %in% paste0("id", c(seq(3, 35, by = 4), 
                                           seq(1, 33, by = 4)))) { ## all dept
        fc_states_sims <- base_sims %>%
          dplyr::mutate(pop_nv = S + E + I + A + R, 
                        pop_1 = S1 + E1 + I1 + A1 + R1, 
                        pop_2 = S2 + E2 + I2 + A2 + R2, 
                        pop_3 = S3 + E3 + I3 + A3 + R3, 
                        pop_4 = S4 + E4 + I4 + A4 + R4,
                        pop_5 = S5 + E5 + I5 + A5 + R5, 
                        pop_6 = S6 + E6 + I6 + A6 + R6,
                        pop_7 = S7 + E7 + I7 + A7 + R7, 
                        pop_8 = S8 + E8 + I8 + A8 + R8,
                        pop_9 = S9 + E9 + I9 + A9 + R9, 
                        pop_10 = S10 + E10 + I10 + A10 + R10,
                        pop = pop_nv + pop_1 + pop_2 + pop_3 + pop_4 + 
                              pop_5 + pop_6 + pop_7 + pop_8 + pop_9 + pop_10, 
                        loglik = start$loglik_f) %>%
          #dplyr::rename(week = time) %>%
          dplyr::rename(sim = .id) %>%
          dplyr::select(sim, week, cases, 
                        S, E, I, A, R, pop_nv, 
                        S1, E1, I1, A1, R1, pop_1, 
                        S2, E2, I2, A2, R2, pop_2,
                        S3, E3, I3, A3, R3, pop_3, 
                        S4, E4, I4, A4, R4, pop_4, 
                        S5, E5, I5, A5, R5, pop_5, 
                        S6, E6, I6, A6, R6, pop_6, 
                        S7, E7, I7, A7, R7, pop_7, 
                        S8, E8, I8, A8, R8, pop_8, 
                        S9, E9, I9, A9, R9, pop_9, 
                        S10, E10, I10, A10, R10, pop_10,  
                        pop, incid, incidU, incidV, asymV, newV, 
                        loglik, foival, Str0, Sout, Sin)
      }
      file_name <- paste0(out_dir_fc, "fcstates_stoch_", fcsettings$str_if2, 
                          "_parid", start$parid, "_seed", fcsettings$seed, 
                          ".csv")
      write_csv(fc_states_sims, file = file_name, ".csv")
    }  %>% 
    group_by(week) %>%
    summarise(S_med = median(S), E_med = median(E), I_med = median(I), 
              A_med = median(A), R_med = median(R), 
              incid_med = median(incid), pop_med = median(pop), 
              cases_med = median(cases),
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
              cases_hi = quantile(cases, probs = c(.975)))
})  -> fc_stoch

plot(fc_stoch$week, fc_stoch$cases_med)

####FORECASTING VISUALIZATIONS##################################################
##----data import and setup-----------------------------------------------------
fig_dir <- paste0("figures/final_", scode, "/")
dir.create(fig_dir, showWarnings = FALSE)

fn_epi <- list.files(out_dir_epi, "fitstates_if_")
fn_epi_full <- paste0(out_dir_epi, fn_epi)
fn_end <- list.files(out_dir_end, "fitstates_if_")
fn_end_full <- paste0(out_dir_end, fn_end)
fn_fc_elim <- list.files(out_dir_fc, "fcstates_stoch_")
fn_fc_elim_full <- paste0(out_dir_fc, fn_fc_elim)
fn_fc_plt <- list.files(out_dir_fc, "stochfcs_")
fn_fc_plt_full <- paste0(out_dir_fc, fn_fc_plt)

haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= episettings$nweeks) 
haiti.dat.end <- haiti.dat %>% dplyr::filter(week >= episettings$nweeks) %>%
  dplyr::mutate(week_end = seq_along(week)-1)

covartab.all <- make.covartab(0, nrow(haiti.dat) + 1, byt = 1, degree = 6, 
                              nbasis = 6, per = 52.14)
covartab <- covartab.all[(episettings$nweeks + 1):nrow(covartab.all), ] %>%
  dplyr::mutate(time_end = seq_along(time) - 1)
num_betas <- ncol(covartab.all) - 1

##----data preparation for plotting---------------------------------------------
if(process_modules){
  fit_epidat <- process_fepi_plot(fn_epi_full, haiti.dat.epi) %>%
    dplyr::mutate(period = "epidemic")
  fit_enddat <- process_fend_plot(fn_end_full, haiti.dat.end) %>% 
    dplyr::filter(week > episettings$nweeks) %>%
    dplyr::mutate(period = "endemic")
  fitdat <- bind_rows(fit_epidat, fit_enddat)
  
  fcdat <- process_ffc_plot(fn_fc_plt_full) %>%
    dplyr::filter(week > max(haiti.dat$week)) %>%
    dplyr::mutate(period = "forecast")
  ffcdat <- bind_rows(fitdat, fcdat)
  
  date.start <- get.mspp.dept.data() %>% 
    dplyr::filter(date_sat == min(date_sat)) %>% 
    dplyr::distinct(date_sat)
  datedat <- data.frame(weekdate = seq(date.start$date_sat,
                                       length.out = nrow(ffcdat), by = 7), 
                        week = seq_len(nrow(ffcdat)))
  
  plotdat_ffc <- full_join(ffcdat, datedat, by = c("week"))
  plotdat_f <- plotdat_ffc %>% dplyr::filter(week <= max(haiti.dat$week))
  plotdat_fc <- plotdat_ffc %>% dplyr::filter(week > max(haiti.dat$week))
}

##----plot clean fits and forecasts---------------------------------------------
if(plot_modules){
  # plot observed fit
  plt_of <- plot_of(plotdat_f)
  ggsave(paste0(fig_dir, "o_ffull_if_", episettings$str_if, "_", 
                endsettings$str_if, ".png"), 
         plt_of, width = 6, height = 4, units = "in")
  print(plt_of)
  
  # plot observed fit without ribbon
  plt_of_nr <- plot_of_nr(plotdat_f)
  ggsave(paste0(fig_dir, "o_ffull_if_", episettings$str_if, "_", 
                endsettings$str_if, "noribbon.png"), 
         plt_of_nr, width = 6, height = 4, units = "in")
  print(plt_of_nr)
  
  # plot observed with zoom
  plt_ofzm <- plot_of(plotdat_f) +
    scale_x_date(limits = c(max(plotdat_f$weekdate) - 300, 
                            max(plotdat_f$weekdate))) +
    scale_y_continuous(limits = c(0, 2000))
  ggsave(paste0(fig_dir, "o_ffullzm_if_", episettings$str_if, "_", 
                endsettings$str_if, ".png"), 
         plt_ofzm, width = 6, height = 4, units = "in")
  print(plt_ofzm)
  
  # plot observed fit with zoom without ribbon
  plt_ofzm_nr <- plot_of_nr(plotdat_f) +
    scale_x_date(limits = c(max(plotdat_f$weekdate) - 300, 
                            max(plotdat_f$weekdate))) +
    scale_y_continuous(limits = c(0, 2000))
  ggsave(paste0(fig_dir, "o_ffullzm_if_", episettings$str_if, "_", 
                endsettings$str_if, "noribbon.png"), 
         plt_ofzm_nr, width = 6, height = 4, units = "in")
  print(plt_ofzm_nr)
  
  # plot observed forecast
  plt_ofc <- plot_ofc(plotdat_fc)
  ggsave(paste0(fig_dir, "o_fcfull_if_", fcsettings$str_if, "_", 
                episettings$str_if, "_", endsettings$str_if, ".png"), 
         plt_ofc, width = 6, height = 4, units = "in")
  print(plt_ofc)
  
  # plot observed forecast without ribbon
  plt_ofc_nr <- plot_ofc_nr(plotdat_fc)
  ggsave(paste0(fig_dir, "o_fcfull_if_", fcsettings$str_if, "_", 
                episettings$str_if, "_", endsettings$str_if, "noribbon.png"), 
         plt_ofc_nr, width = 6, height = 4, units = "in")
  print(plt_ofc_nr)
  
  # plot observed fit and forecast
  plt_offc <- plot_offc_full(plotdat_ffc)
  ggsave(paste0(fig_dir, "o_ffcfull_if_", fcsettings$str_if, "_", 
                episettings$str_if, "_", endsettings$str_if, ".png"), 
         plt_offc, width = 6, height = 4, units = "in")
  print(plt_offc)
  
  # plot observed fit and forecast without ribbon
  plt_offc_nr <- plot_offc_full_nr(plotdat_ffc)
  ggsave(paste0(fig_dir, "o_ffcfull_if_", fcsettings$str_if, "_", 
                episettings$str_if, "_", endsettings$str_if, "noribbon.png"), 
         plt_offc_nr, width = 6, height = 4, units = "in")
  print(plt_offc_nr)
  
  # plot true forecast
  plt_tfc <- plot_tfc(plotdat_fc)
  ggsave(paste0(fig_dir, "t_fcfull_if_", fcsettings$str_if, "_", 
                episettings$str_if, "_", endsettings$str_if, ".png"), 
         plt_tfc, width = 6, height = 4, units = "in")
  print(plt_tfc)
  
  # plot true forecast without ribbon
  plt_tfc_nr <- plot_tfc_nr(plotdat_fc)
  ggsave(paste0(fig_dir, "t_fcfull_if_", fcsettings$str_if, "_", 
                episettings$str_if, "_", endsettings$str_if, "noribbon.png"), 
         plt_tfc_nr, width = 6, height = 4, units = "in")
  print(plt_tfc_nr)
  
  # plot true fit and forecast
  plt_tffc <- plot_tffc_full(plotdat_ffc)
  ggsave(paste0(fig_dir, "t_ffcfull_if_", fcsettings$str_if, "_", 
                episettings$str_if, "_", endsettings$str_if, ".png"), 
         plt_tffc, width = 6, height = 4, units = "in")
  print(plt_tffc)
  
  # plot true fit and forecast without ribbon
  plt_tffc_nr <- plot_tffc_full_nr(plotdat_ffc)
  ggsave(paste0(fig_dir, "t_ffcfull_if_", fcsettings$str_if, "_", 
                episettings$str_if, "_", endsettings$str_if, "noribbon.png"), 
         plt_tffc_nr, width = 6, height = 4, units = "in")
  print(plt_tffc_nr)
}

####MODEL FIT VISUALIZATIONS####################################################
##----plot cumulative vaccination across scenarios------------------------------
if(cumvac_modules) {
  cumvac_df <- process_cumvac_plot(datedat)
  plt_cumvac <- plot_cumvac(cumvac_df)
  ggsave(paste0(fig_dir, "cumvacc_deployment.png"), 
         plt_cumvac, width = 6, height = 4, units = "in")
  plt_cumvac_frac <- plot_cumvac_frac(cumvac_df, plotdat_fc)
  ggsave(paste0(fig_dir, "cumvacc_deployment_frac.png"), 
         plt_cumvac_frac, width = 6, height = 4, units = "in")
  print(plt_cumvac)
  print(plt_cumvac_frac)
}

##----data preparation for elimination calculations-----------------------------
if(!file.exists(paste0(out_dir_fc, "fcelim_if_", fcsettings$str_if, ".rds"))) {
  fcelim <- process_fc_elim(fn_fc_elim_full, scen_code)
  saveRDS(fcelim, paste0(out_dir_fc, "fcelim_if_", fcsettings$str_if, ".rds"))
} else {
  fcelim <- readRDS(paste0(out_dir_fc, "fcelim_if_", fcsettings$str_if, ".rds"))
}

##----print elimination---------------------------------------------------------
pElim_dat <- find_probElim(fcelim)
print(paste("probability of elimination", fcsettings$str_if))
print(pElim_dat)

tElim_dat <- find_timeToElim(fcelim)
print(paste("time to elimination (weeks)"))
print(tElim_dat)

sink(paste0(out_dir_fc, "fcelimprob.txt"))
print.table(pElim_dat)
print(paste("*********"))
print.table(tElim_dat)
sink()
