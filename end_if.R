#==============================================================================#
# Model 1: Endemic Iterated Filtering                                          #
#          Parameter estimation for endemic period (Apr. 2015 - Jan. 2019)     #
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

####ENDEMIC MODEL FITTING#######################################################
#### settings and data
episettings <- fix_epi_settings(nprofs = 3, nparticles = numparts, nmif = 150)
endsettings <- fix_end_settings(episettings)
fig_dir <- paste0("figures/", endsettings$mcode, "/")
out_dir_epi <- paste0("GeneratedData/", episettings$mcode, "/")
out_dir_end <- paste0("GeneratedData/", endsettings$mcode, "/")
out_dir <- paste0("GeneratedData/", endsettings$mcode, "/")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir_end, showWarnings = FALSE)

starts <- generate_end_params(out_dir_epi = out_dir_epi, out_dir_end = out_dir_end, 
                              episettings = episettings, endsettings = endsettings)
# pop.haiti <- get_haiti_pop()
# 
# haiti.dat <- get_mspp_agg_data() 
# haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= episettings$nweeks) 
# haiti.dat.end <- haiti.dat %>% dplyr::filter(week > episettings$nweeks) %>%
#   dplyr::mutate(week_end = seq_along(week))

# covartab <- make_covartab(0, nrow(haiti.dat) + 1, byt = 1, ## make table of covariates for seasonality
#                           degree = 6, nbasis = 6, per = 52.14)

# haiti.mod <- mod_build(dat = haiti.dat, covar = covartab) ## build model

#### perform mif2 on endemic period
prof_if <- mod_fit(endsettings, out_dir, 1, starts, haiti.mod)

####ENDEMIC MODEL FITTING CHECKS################################################
#### look at just one set of parameter estimates
ismp <- sample(1:nrow(prof_if), size=1)

#### record parameter estimates as .csv
val_ests <- rec_res(prof_if)
write_csv(val_ests, paste0(out_dir, "end_par_ests", endsettings$str_if2,
                           "_seed", endsettings$seed, ".csv"))

#### check beta and rough R0 
betas <- prof_if[ismp,] %>% dplyr::select(contains("beta")) %>% as.matrix
bspline <- covartab %>% dplyr::select(contains("seas")) %>% as.matrix %>% t
mybeta <- betas %*% bspline %>% as.vector
plot(mybeta, type = "l")

#### adjustment for mixing ???
I0s <- prof_if[ismp,] %>% dplyr::select(I_0) %>% unlist %>% median
nus <- prof_if[ismp,] %>% dplyr::select(nu) %>% unlist %>% median
print(paste("R0 ranges (init, min, max):", (I0s^nus)/I0s*mybeta[1]*2/7, 
            (I0s^nus)/I0s*min(mybeta)*2/7, (I0s^nus)/I0s*max(mybeta)*2/7))

#### explore parameter estimates and likelihoods
pltcols <-  c("rho", "tau", "beta1", "nu", "loglik")
all_if <- prof_if %>% dplyr::select(pltcols)
plt_all_if <- ggpairs(all_if)
print(plt_all_if)

#### get one set of estimates and simulate data for epidemic period
allpars <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6", ## all parameter names
             "gamma","sigma","theta0","alpha","mu","delta","nu",
             "S_0","E_0","I_0","A_0","R_0", "pop_0")
params <- parmat(unlist(prof_if[ismp, colnames(prof_if) %in% allpars]), nrep=1) ## get values
end_mod <- window(haiti.mod, start = 233, end = 430) ## subset times from full model
timezero(end_mod) <- 232

simdat <- sim_data(end_mod, params, 25) %>% ## simulate and join with case data
  dplyr::select(week, cases_mean) %>%
  dplyr::rename(cases = cases_mean) %>%
  dplyr::left_join(haiti.dat %>% dplyr::rename(orig=cases), by = c("week"))

#### plot simulated data 
fitplt <- plot_sim_dat(simdat, period = 1)
print(fitplt)

#### plot simulated data from all sets of parameter estimates
fn_if <- list.files(out_dir, "fitstates_if_")
fn_if_full <- paste0(out_dir, fn_if)
fit_if_summ <- process_fitstates(fn_if_full, haiti.dat.end)

plt_if_summ_wribbon <- plot_fit_wribbon(fit_if_summ)
print(plt_if_summ_wribbon)

plt_if_summ_noribbon <- plot_fit_noribbon(fit_if_summ)
print(plt_if_summ_noribbon)

# plt_if_summ_zm_wribbon <- plot_fit_wribbon(fit_if_summ) +
#   scale_x_continuous(limits = c(300, 380)) +
#   scale_y_continuous(limits = c(0, 3000))
# print(plt_if_summ_zm_wribbon)
# 
# plt_if_summ_zm_noribbon <- plot_fit_noribbon(fit_if_summ) +
#   scale_x_continuous(limits = c(300, 380)) +
#   scale_y_continuous(limits = c(0, 3000))
# print(plt_if_summ_zm_noribbon)

#### explore full epidemic/endemic iterated filtering fits
fn_epi <- list.files(out_dir_epi, "fitstates_if_")
fn_epi_full <- paste0(out_dir_epi, fn_epi)
fn_end <- list.files(out_dir_end, "fitstates_if_")
fn_end_full <- paste0(out_dir_end, fn_if)

fit_epidat <- process_fitstates(fn_epi_full, haiti.dat.epi) %>%
  dplyr::mutate(period = "epidemic")
fit_enddat <- process_fitstates(fn_end_full, haiti.dat.end) %>%
  dplyr::filter(week > episettings$nweeks) %>%
  dplyr::mutate(period = "endemic")
fitdat <- bind_rows(fit_epidat, fit_enddat) %>%
  dplyr::select(-week_end)

plt_fullfit_wribbon <- plot_fit_full_wribbon(fitdat)
print(plt_fullfit_wribbon)

plt_fullfit_noribbon <- plot_fit_full_noribbon(fitdat)
print(plt_fullfit_noribbon)

# plt_fullfit_zm_wribbon <- plt_fullfit_wribbon +
#   scale_x_continuous(limits = c(max(fitdat$week)-40, max(fitdat$week))) +
#   scale_y_continuous(limits = c(0, 2000))
# print(plt_fullfit_zm_wribbon)
# 
# plt_fullfit_zm_noribbon <- plt_fullfit_noribbon +
#   scale_x_continuous(limits = c(max(fitdat$week)-40, max(fitdat$week))) +
#   scale_y_continuous(limits = c(0, 2000))
# print(plt_fullfit_zm_noribbon)

#### try simulating far into the future to test forecasting
fcsettings <- fix_fc_vac_settings("id0", nsims = 25)

fig_dir <- paste0("figures/", fcsettings$mcode, "_", fcsettings$scode, "/")
out_dir_fc <- paste0("GeneratedData/", fcsettings$mcode, "_", fcsettings$scode, "/")
out_dir_end <- paste0("GeneratedData/", endsettings$mcode, "/")
out_dir_epi <- paste0("GeneratedData/", episettings$mcode, "/")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir_fc, showWarnings = FALSE)
colvec <- c("endemic" = "#8d006a", "forecast" = "#a890ac")

## predict into future with my params
fc_test_mod <- haiti.mod
time(fc_test_mod, t0 = TRUE) <- seq(232, 1003)
fc_simdat <- sim_data(fc_test_mod, params, num_sims = 25) %>%
  rename(est_lo = cases_lo) %>%
  rename(est_hi = cases_hi) %>%
  rename(est_med = cases_med)
fc_simdat$orig_cases <- NA
fc_simdat$week_end <- NA
fc_simdat <- fc_simdat %>%
  mutate(period = ifelse(week < 433, "endemic", "forecast"))
fc_simdat <- fc_simdat %>%
  select(week, est_med, est_lo, est_hi, orig_cases, week_end, period)
fc_simdat$date <- lubridate::ymd("2010-10-23") + lubridate::weeks(fc_simdat$week)

plt1 <- ggplot(fc_simdat, aes(x = date)) +
  geom_line(aes(y = est_med, colour = period), size = 1) +
  geom_ribbon(aes(ymin = est_lo, ymax = est_hi, fill = period), alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_labels = "%Y %b", date_breaks = "1 year") +
  scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))

## predict into future with supplement params
allpars <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6",
             "gamma","sigma","theta0","alpha","mu","delta","nu",
             "S_0","E_0","I_0","A_0","R_0", "pop_0")
pars <- c(0.34, 4.1, 5.4, 3.3, 4.3, 3.5, 5.0, 2.7,
          7/2, 7/1.4, 0, 7/2920, (((1+22.6/1000)^(1/52.14))-1), (((1+7.5/1000)^(1/52.14))-1), 0.96, 
          1 - 248 / 11658690 - 324 / 11658690 - 417522 / 11658690,
          248 / 11658690, 324 / 11658690, 0.0, 1417522 / 11658690, 11658690)
names(pars) <- allpars

fc_simdat <- sim_data(fc_test_mod, pars, num_sims = 25) %>%
  rename(est_lo = cases_lo) %>%
  rename(est_hi = cases_hi) %>%
  rename(est_med = cases_med)
fc_simdat$orig_cases <- NA
fc_simdat$week_end <- NA
fc_simdat <- fc_simdat %>%
  mutate(period = ifelse(week < 433, "endemic", "forecast"))
fc_simdat <- fc_simdat %>%
  select(week, est_med, est_lo, est_hi, orig_cases, week_end, period)
fc_simdat$date <-lubridate::ymd("2010-10-23") + lubridate::weeks(fc_simdat$week)

#### plot simulated data 
plt2 <- ggplot(fc_simdat, aes(x = date)) +
  geom_line(aes(y = est_med, colour = period), size = 1) +
  geom_ribbon(aes(ymin = est_lo, ymax = est_hi, fill = period), alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_labels = "%Y %b", date_breaks = "1 year") +
  scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))

plt1
plt2
