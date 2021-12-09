#==============================================================================#
# Model 1: Epidemic Iterated Filtering                                         #
#          Parameter estimation for epidemic period (Oct. 2010 - Mar. 2015)    #
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
num_iters <- 150 # number of iterations for IF2 (originally 100)
num_parts <- 5000 # number of particles for IF2 (originally 100)
num_sims <- 15 # number of simulations for median calculations (originally 25)

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

####EPIDEMIC MODEL FITTING######################################################
#### settings and data
episettings <- fix_epi_settings(nprofs = 3, nparticles = num_parts, nmif = num_iters)
endsettings <- fix_end_settings(episettings)
fcsettings <- fix_fc_vac_settings("id0", nsims = 25)
fig_dir <- paste0("figures/", episettings$mcode, "/")
out_dir <- paste0("GeneratedData/", episettings$mcode, "/")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)

starts <- generate_epi_params(out_dir = out_dir, departement = 'Artibonite') ## generate starting values for parameter estimation
pop.haiti <- get_haiti_pop(departement = 'Artibonite') ## get haiti population (country-wide)
# haiti.dat <- get_mspp_agg_data() ## get case data (country-wide)
haiti.dat <- haiti1_data()
haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= episettings$nweeks) ## get case data (epidemic period)
haiti.dat.end <- haiti.dat %>% dplyr::filter(week > episettings$nweeks) %>% ## get case data (endemic period)
  dplyr::mutate(week_end = seq_along(week))
# covartab <- make_covartab(0, nrow(haiti.dat)+fcsettings$horizon+1,
#                               byt=1, degree=6, nbasis=6, per=52.14)
# haiti.mod <- mod_build(dat = haiti.dat, covar = covartab) ## build model
# haiti.mod <- haiti1()
haiti2 <- haiti1_dep(departement = "Artibonite", vacscen = "id0")

#### perform mif2 on epidemic period
prof_if <- mod_fit(episettings, out_dir, 0, starts, haiti2)

####EPIDEMIC MODEL FITTING CHECKS###############################################
#### look at just one set of parameter estimates
ismp <- sample(1:nrow(prof_if), size=1)

#### record parameter estimates as .csv
val_ests <- rec_res(prof_if)
write_csv(val_ests, paste0(out_dir, "epi_par_ests", episettings$str_if2,
                           "_seed", episettings$seed, ".csv"))

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
pltcols <-  c("rho", "tau", "beta1", "nu", "sig_sq", "loglik")
all_if <- prof_if %>% dplyr::select(pltcols)
plt_all_if <- ggpairs(all_if)
print(plt_all_if)

#### get one set of estimates and simulate data for epidemic period
allpars <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6", ## all parameter names
             "gamma","sigma","theta0","alpha","mu","delta","nu", "sig_sq",
             "S_0","E_0","I_0","A_0","R_0", "pop_0")
params <- parmat(unlist(prof_if[ismp, colnames(prof_if) %in% allpars]), nrep=1) ## get values
epi_mod <- window(haiti.mod, start = 0, end = 232) ## subset times from full model

simdat <- sim_data(epi_mod, params, 25) %>% ## simulate and join with case data
  dplyr::select(week, cases_mean) %>%
  dplyr::rename(cases = cases_mean) %>%
  dplyr::left_join(haiti.dat %>% dplyr::rename(orig=cases), by = c("week"))

#### plot simulated data
fitplt <- plot_sim_dat(simdat, period = 0, prof_if)
print(fitplt)

#### plot simulated data from all sets of parameter estimates
fn_if <- list.files(out_dir, "fitstates_if_")
fn_if_full <- paste0(out_dir, fn_if)
fit_if_summ <- process_fitstates(fn_if_full, haiti.dat.epi)

plt_if_summ_wribbon <- plot_fit_wribbon(fit_if_summ)
print(plt_if_summ_wribbon)

plt_if_summ_noribbon <- plot_fit_noribbon(fit_if_summ)
print(plt_if_summ_noribbon)

plt_if_summ_zm_wribbon <- plot_fit_wribbon(fit_if_summ) +
  scale_x_continuous(limits = c(episettings$nweeks-80,episettings$nweeks)) +
  scale_y_continuous(limits = c(0, 6000))
print(plt_if_summ_zm_wribbon)

plt_if_summ_zm_noribbon <- plot_fit_noribbon(fit_if_summ) +
  scale_x_continuous(limits = c(episettings$nweeks-80,episettings$nweeks)) +
  scale_y_continuous(limits = c(0, 6000))
print(plt_if_summ_zm_noribbon)
