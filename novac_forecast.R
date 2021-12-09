#==============================================================================#
# Model 1: No Vaccination Forecasting                                          #
#          Predicting resurgence and elimination without a vaccination         #
#          campaign                                                            #
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

####NO VAC MODEL FITTING########################################################
episettings <- fix_epi_settings(nprofs = 3, nparticles = 1000, nmif = 150)
endsettings <- fix_end_settings(episettings)
fcsettings <- fix_fc_vac_settings("id0", nsims = num_sims)

fig_dir <- paste0("figures/", fcsettings$mcode, "_", fcsettings$scode, "/")
out_dir_fc <- paste0("GeneratedData/", fcsettings$mcode, "_", fcsettings$scode, "/")
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

covartab.all <- make_covartab(0, nrow(haiti.dat)+fcsettings$horizon+1, 
                              byt=1, degree=6, nbasis=6, per=52.14) %>%
  dplyr::mutate(time_endfc = ifelse(time < episettings$nweeks, NA, time)) %>%
  dplyr::mutate(time_fc = ifelse(time < max(haiti.dat$week), NA, time)) %>%
  dplyr::mutate(time_end = ifelse(time <= max(haiti.dat$week), time, NA))
covartab.all2 <- covartab.all %>%
  dplyr::select(time, contains("seas"))
num_betas <- ncol(covartab.all2) - 1

##----make models-------------------------------------------------------------##
full_mod <- mod_fc_build(dat = haiti.dat, covar = covartab.all2,
                         depts = depts)
haiti.mod.fc <- full_mod

##----perform stochastic forecasts--------------------------------------------##
bake(file=paste0(out_dir_fc, "/stochfcs_", fcsettings$str_if, ".rds"), 
     seed = fcsettings$seed, {
  foreach(start = iter(starts, by = 'row'),
          .inorder = FALSE, .combine = rbind,
          .packages = c("pomp", "magrittr"),
          .errorhandling = c('remove'),
          .noexport = c(),
          .verbose = TRUE) %dopar%
    {
      tic <- Sys.time()
      print(paste("forecast stoch params novac", start$parid))
      
      M2 <- haiti.mod.fc
      time(M2) <- seq((tail(haiti.dat.end$week, n = 1) + 1), tail(covartab.all2$time, n = 1))
      timezero(M2) <- tail(haiti.dat.end$week, n = 1)
      allpars <- c("rho", "tau", "beta1", "beta2", "beta3", "beta4", 
                   "beta5", "beta6", "gamma", "sigma", "theta0", "alpha",
                   "mu", "delta", "nu", "S_0", "E_0", "I_0", "A_0",
                   "R_0", "incid_0", "pop_0")
      coef(M2) <- unlist(start[which(names(start) %in% allpars)])
      fc_states_sims <- sim_data(M2, num_sim = 25)
      write_csv(fc_states_sims, paste0(out_dir_fc, "fcstates_stoch_", 
                                       fcsettings$str_if2, "_parid", 
                                       start$parid, "_seed", 
                                       fcsettings$seed, ".csv"))
    } %>% group_by(week) %>% ## medians of medians
         summarise(S_med = median(S_med), E_med = median(E_med), I_med = median(I_med), 
                   A_med = median(A_med), R_med = median(R_med), 
                   incid_med = median(incid_med), pop_med = median(pop_med), 
                   cases_med = median(cases_med),
                   S_lo = quantile(S_med, probs = c(.025)), 
                   E_lo = quantile(E_med, probs = c(.025)), 
                   I_lo = quantile(I_med, probs = c(.025)), 
                   A_lo = quantile(A_med, probs = c(.025)), 
                   R_lo = quantile(R_med, probs = c(.025)), 
                   incid_lo = quantile(incid_med, probs = c(.025)), 
                   pop_lo = quantile(pop_med, probs = c(.025)), 
                   cases_lo = quantile(cases_med, probs = c(.025)),
                   S_hi = quantile(S_med, probs = c(.975)), 
                   E_hi = quantile(E_med, probs = c(.975)), 
                   I_hi = quantile(I_med, probs = c(.975)),
                   A_hi = quantile(A_med, probs = c(.975)), 
                   R_hi = quantile(R_med, probs = c(.975)), 
                   incid_hi = quantile(incid_med, probs = c(.975)), 
                   pop_hi = quantile(pop_med, probs = c(.975)), 
                   cases_hi = quantile(cases_med, probs = c(.975)))
})  -> fc_stoch

####FORECASTING VISUALIZATIONS##################################################
##----data import and setup---------------------------------------------------##
fn_epi <- list.files(out_dir_epi, "fitstates_if_")
fn_epi_full <- paste0(out_dir_epi, fn_epi)
fn_end <- list.files(out_dir_end, "fitstates_if_")
fn_end_full <- paste0(out_dir_end, fn_end)
fn_fc_elim <- list.files(out_dir_fc, "fcstates_stoch_")
fn_fc_elim_full <- paste0(out_dir_fc, fn_fc_elim)
fn_fc_plt <- list.files(out_dir_fc, "stochfcs_")
fn_fc_plt_full <- paste0(out_dir_fc, fn_fc_plt)

##----data preparation for plotting-------------------------------------------##
if(plot_modules) {
  fit_epidat <- process_fc_plot(fn_epi_full, haiti.dat.epi) %>%
    dplyr::mutate(period = "epidemic")
  fit_enddat <- process_fc_plot(fn_end_full, haiti.dat.end) %>%
    dplyr::mutate(period = "endemic")
  fitdat <- bind_rows(fit_epidat, fit_enddat)
  
  fcdat <- process_ffc_plot(fn_fc_plt_full) %>%
    dplyr::filter(week > max(haiti.dat$week)) %>%
    dplyr::mutate(period = "forecast")
  ffcdat <- bind_rows(fitdat, fcdat)
  
  date.start <- get_mspp_dept_data() %>% 
    dplyr::filter(date_sat == min(date_sat)) %>% 
    dplyr::distinct(date_sat)
  datedat <- data.frame(weekdate = seq(date.start$date_sat, 
                                       length.out = nrow(ffcdat), by = 7), 
                        week = seq_len(nrow(ffcdat)))
  
  plotdat_ffc <- full_join(ffcdat, datedat, by = c("week")) ## fit and forecast
  plotdat_f <- plotdat_ffc %>% dplyr::filter(week <= max(haiti.dat$week)) ## fit
  plotdat_fc <- plotdat_ffc %>% dplyr::filter(week > max(haiti.dat$week)) ## forecast
}

##----plot clean fits and forecasts-------------------------------------------##
if(plot_modules){
  # plot fit with case data
  plt_of <- plot_fc_case(plotdat_f)
  print(plt_of)
  
  # plot forecast - reported cases
  plt_ofc <- plot_fc_case(plotdat_fc)
  print(plt_ofc)
  
  # plot fit and forecast - reported cases
  plt_offc <- plot_fc_case(plotdat_ffc)
  print(plt_offc)
  
  # plot forecast - incidence
  plt_tfc <- plot_fc_incid(plotdat_fc)
  print(plt_tfc)
  
  # plot fit and forecast - incidence
  plt_tffc <- plot_fc_incid(plotdat_ffc)
  print(plt_tffc)
}

####MODEL FIT VISUALIZATIONS####################################################
##----data preparation for elimination calculations---------------------------##
if (!file.exists(paste0(out_dir_fc, "fcelim_if_", fcsettings$str_if, ".rds"))) {
  fcelim <- process_fc_elim(fcsettings, fn_fc_elim_full)
  saveRDS(fcelim, paste0(out_dir_fc, "fcelim_if_", fcsettings$str_if, ".rds"))
} else {
  fcelim <- readRDS(paste0(out_dir_fc, "fcelim_if_", fcsettings$str_if, ".rds"))
}

##----print elimination-------------------------------------------------------##
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
