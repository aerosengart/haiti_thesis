##=============================================================================##
## Reproduction: Thesis Code                                                   ##
##=============================================================================##
## load packages/functions
library("ggplot2")
library("pomp")
library("dplyr")
library("purrr")
library("tidyverse")
library("magrittr")
library("GGally")
library("gridExtra")
library("cowplot")
library("doParallel")
source("mod_tools.R")
source("mod_vis.R")
source("forecast_tools.R")
source("mod_build.R")
source("haiti1_vacscen.R")
source("haiti1.R")
source("haiti1_joint.R")
source("haiti1_agg_data.R")
source("haiti1_data.R")
source("haiti1_covars.R")

plot_sim_dat <- function(simdat) {
  plt <- ggplot(simdat, aes(x = week)) +
    geom_line(aes(y = cases), colour = "blue", alpha = 0.5) +
    geom_point(aes(y = orig), size = 1) +
    theme_bw() +
    xlim(0, 430) +
    scale_y_continuous("cases", limits = c(0, 30000), oob = scales::squish)
  return(plt)
}

haiti_dat <- haiti1_agg_data()

################################################################################
## ARMA Benchmarks                                                            ##
################################################################################
## epidemic period
epi_haiti_dat <- haiti_dat %>% filter(week <= 232)
epi_log_full_cases <- log(epi_haiti_dat$cases + 1)
epi_log_sum <- sum(epi_log_full_cases)
epi_arima_21 <- arima(epi_log_full_cases, order = c(2, 0, 1))
epi_arima_21$loglik - epi_log_sum

## endemic period
end_haiti_dat <- haiti_dat %>% filter(week > 232)
end_log_full_cases <- log(end_haiti_dat$cases + 1)
end_log_sum <- sum(end_log_full_cases, na.rm = TRUE)
end_arima_21 <- arima(end_log_full_cases, order = c(2, 0, 1))
end_arima_21$loglik - end_log_sum 

## full period
log_full_cases <- log(haiti_dat$cases + 1)
log_sum <- sum(log_full_cases, na.rm = TRUE)
arima_21 <- arima(log_full_cases, order = c(2, 0, 1))
arima_21$loglik - log_sum 

################################################################################
## Reproduction - Figure A1                                                   ##
################################################################################
load("data/orig_mod/epi/epidemic_fit.rda")
epi_pars <- prof_if
parsepi_3000 <- epi_pars %>%
  dplyr::filter(nu > 0.9) %>%
  dplyr::filter(beta1 < 100 & beta2 < 100 & beta3 < 100 & beta4 < 100 & beta5 < 100 & beta6 < 100) 
parids_epi <- parsepi_3000$parid %>% as.data.frame()
load("data/orig_mod/end/orig_end_mif.rda")
end_pars <- prof_if
parsend_3000 <- end_pars %>%
  dplyr::filter(loglik > -3000)
parids_end <- parsend_3000$parid %>% as.data.frame()
parids <- dplyr::inner_join(parids_end, parids_epi)
parids <- parids$.
epi_recs_fullpath <- c(paste0("data/orig_mod/epi/fitstates_if_epi232_nprofs3_nmif100_nparticles100_parid",
                       parids, "_seed20190415.csv"))
end_recs_fullpath <- c(paste0("data/orig_mod/end/fitstates_if_parid",
                       parids, ".csv"))
simdat_epi <- process_fitstates(epi_recs_fullpath, haiti_dat_epi)
simdat_end <- process_fitstates(end_recs_fullpath, haiti_dat_end) %>%
  dplyr::filter(week > 232) %>%
  dplyr::select(-week_end)
simdat_all <- rbind(simdat_epi, simdat_end)
simdat_all$period <- ifelse(simdat_all$week < 233, "epidemic", "endemic")
simdat_all_log <- simdat_all %>%
  dplyr::mutate(orig_cases = log(orig_cases  + 1),
                est_med = log(est_med + 1),
                est_lo = log(est_lo + 1),
                est_hi = log(est_hi + 1))
reprod_plot <- plot_fit_full_wribbon(simdat_all) + ggtitle("Reproduction Plot: Natural Scale")
reprod_plot_log <- plot_fit_full_wribbon(simdat_all_log) + ggtitle("Reproduction Plot: Log Scale") + ylab("log(Reported Cases + 1)")
ggsave(reprod_plot_log, file = "reprod_plot_log.png", width = 8, height = 6)
inset_dat <- simdat_all %>%
  dplyr::filter(week > 115)
inset <- plot_fit_full_wribbon(inset_dat) + theme(legend.position = "none", 
                                                  axis.title.x=element_blank(), 
                                                  axis.title.y=element_blank(), 
                                                  axis.text.x = element_text(size = 6),
                                                  axis.text.y = element_text(size = 6))
reprod_inset <-
  ggdraw() +
  draw_plot(reprod_plot) +
  draw_plot(inset, x = .4, y = .5, width = .55, height = .4)
ggsave(reprod_inset, file = "reprod_inset.png", width = 8, height = 6)


## Reproduction Likelihood
best_pars_epi_reprod <- parsepi_3000[which.max(parsepi_3000$loglik), ]
best_pars_end_reprod <- parsend_3000[which.max(parsend_3000$loglik), ]
allpars <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6", 
             "gamma","sigma","theta0","alpha","mu","delta","nu", 
             "S_0","E_0","I_0","A_0","R_0", "pop_0")
reprod_epi_params <- best_pars_epi_reprod[colnames(best_pars_epi_reprod) %in% allpars]
reprod_end_params <- best_pars_end_reprod[colnames(best_pars_end_reprod) %in% allpars]
episettings <- fix_epi_settings(nprofs = 3, nparticles = 5000, nmif = 150)
fcsettings <- fix_fc_vac_settings("id0", nsims = 25)
## get case data (epidemic period)
haiti_dat_epi <- haiti_dat %>% dplyr::filter(week <= episettings$nweeks)
## get case data (endemic period)
haiti_dat_end <- haiti_dat %>% dplyr::filter(week > episettings$nweeks) %>%
  dplyr::mutate(week_end = seq_along(week))
## make table of covariates for seasonality
covartab <- make_covartab(0, nrow(haiti_dat) + fcsettings$horizon + 1, byt = 1,
                          degree = 6, nbasis = 6, per = 52.14)
full_mod <- mod_build(dat = haiti_dat, covar = covartab)
epi_mod <- full_mod %>%
  pomp(params = reprod_epi_params) %>%
  window(start = 0, end = episettings$nweeks)
end_mod <- full_mod %>%
  window(start = episettings$nweeks + 1, end = nrow(haiti_dat)) %>%
  pomp(params = reprod_end_params, t0 = episettings$nweeks)
replicate(20, epi_mod %>% pfilter(Np = 5000) %>% logLik()) %>%
  logmeanexp(se=TRUE) -> ll_epi_reprod
replicate(20, end_mod %>% pfilter(Np = 5000) %>% logLik()) %>%
  logmeanexp(se=TRUE) -> ll_end_reprod

################################################################################
## Adjusted Model Simulation Plot - Figure 2                                  ##
################################################################################
load(file = "data/pruned_epi.rda")
epi_par_df <- pruned_epi
best_pars <- epi_par_df[which.max(epi_par_df$loglik), ]
allpars <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6", 
             "gamma","sigma","theta0","alpha","mu","delta","nu", "sig_sq",
             "S_0","E_0","I_0","A_0","R_0", "pop_0")
adj_epi_params <- best_pars[colnames(best_pars) %in% allpars]
epi_mod <- haiti1()
simdat_mle_epi <- sim_data(epi_mod, adj_epi_params, 25)
load(file = "data/pruned_end.rda")
end_par_df <- pruned_end
best_pars <- end_par_df[which.max(end_par_df$loglik), ]
allpars <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6",
             "gamma","sigma","theta0","alpha","mu","delta","nu", "sig_sq",
             "S_0","E_0","I_0","A_0","R_0", "pop_0")
adj_end_params <- best_pars[colnames(best_pars) %in% allpars]
end_mod <- haiti1(period = "endemic")
simdat_mle_end <- sim_data(end_mod, adj_end_params, 25)
simdat_mle <- rbind(simdat_mle_epi, simdat_mle_end)
simdat_mle <- simdat_mle %>%
  dplyr::left_join(haiti_dat, by = c("week")) %>%
  dplyr::rename(est_lo = cases_lo,
                est_hi = cases_hi,
                est_med = cases_med,
                orig_cases = cases)
simdat_mle$period <- ifelse(simdat_mle$week < 233, "epidemic", "endemic")
simdat_mle_log <- simdat_mle %>%
  dplyr::mutate(orig_cases = log(orig_cases  + 1),
                est_med = log(est_med + 1),
                est_lo = log(est_lo + 1),
                est_hi = log(est_hi + 1))
mle_plot_adj <- plot_fit_full_wribbon(simdat_mle) + ggtitle("Adjusted Model Plot: Natural Scale")
mle_plot_log <- plot_fit_full_wribbon(simdat_mle_log) + ggtitle("Adjusted Model Plot: Log Scale") + ylab("log(Reported Cases + 1)")
ggsave(mle_plot_log, file = "mle_plot_adj_log.png", width = 8, height = 6)
inset_dat <- simdat_mle %>%
  dplyr::filter(week > 115)
inset <- plot_fit_full_wribbon(inset_dat) + theme(legend.position = "none", 
                                                  axis.title.x=element_blank(), 
                                                  axis.title.y=element_blank(), 
                                                  axis.text.x = element_text(size = 6),
                                                  axis.text.y = element_text(size = 6))
adj_inset <-
  ggdraw() +
  draw_plot(mle_plot_adj) +
  draw_plot(inset, x = .4, y = .5, width = .55, height = .4)
ggsave(adj_inset, file = "adj_inset.png", width = 8, height = 6)

## Adjusted Likelihood
coef(epi_mod) <- adj_epi_params
replicate(20, epi_mod %>% pfilter(Np = 5000) %>% logLik()) %>%
  logmeanexp(se=TRUE) -> ll_epi_adj
coef(end_mod) <- adj_end_params
replicate(20, end_mod %>% pfilter(Np = 5000) %>% logLik()) %>%
  logmeanexp(se=TRUE) -> ll_end_adj

################################################################################
## Profile Likelihood - Figure 3                                              ##
################################################################################
load("data/profile_likelihood_sigsq_epi.rda")
results %>%
  filter(loglik>max(loglik)-100,loglik.se<1) %>%
  group_by(round(sig_sq,2)) %>%
  filter(rank(-loglik)<3) %>%
  ungroup() %>%
  ggplot(aes(x=sig_sq,y=loglik))+
  geom_point()+
  labs(title = "epidemic") +
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  ) -> epi_sig_sq
load("data/profile_likelihood_sigsq_end.rda")
results_end <- results %>%
  filter(!is.na(loglik))
results_end %>%
  filter(!is.na(loglik)) %>%
  filter(loglik>max(loglik)-100, loglik.se<1) %>%
  group_by(round(sig_sq,2)) %>%
  filter(rank(-loglik)<3) %>%
  ungroup() %>%
  ggplot(aes(x=sig_sq,y=loglik))+
  geom_point()+
  labs(title = "endemic") +
  geom_hline(
    color="red",
    yintercept=max(results_end$loglik)-0.5*qchisq(df=1,p=0.95)
  ) -> end_sig_sq
prof_lik_plot <- arrangeGrob(epi_sig_sq, end_sig_sq, ncol = 2)

################################################################################
## Adjusted, Joint Model Plot - Figure 4                                      ##
################################################################################
load("data/joint_mif_ssqtaurho.rda")
joint_mod <- haiti1_joint()
best_pars <- prof_if[which.max(prof_if$full_loglik), ] %>%
  select(-c(1:2, 29:34)) %>%
  rename(beta1 = beta1_joint,
         beta2 = beta2_joint,
         beta3 = beta3_joint,
         beta4 = beta4_joint,
         beta5 = beta5_joint,
         beta6 = beta6_joint)
simdat_test <- sim_data(joint_mod, best_pars, 50) %>%
  dplyr::left_join(haiti_dat, by = c("week")) %>%
  dplyr::rename(est_lo = cases_lo,
                est_hi = cases_hi,
                est_med = cases_med,
                orig_cases = cases)
simdat_test$period <- ifelse(simdat_test$week < 233, "epidemic", "endemic")
simdat_mle_log <- simdat_test %>%
  dplyr::mutate(orig_cases = log(orig_cases  + 1),
                est_med = log(est_med + 1),
                est_lo = log(est_lo + 1),
                est_hi = log(est_hi + 1))


mle_plot_joint <- plot_fit_full_wribbon(simdat_test) + ggtitle("Adjusted, Joint Model Plot: Natural Scale")
mle_plot_joint_log <- plot_fit_full_wribbon(simdat_mle_log) + ggtitle("Adjusted, Joint Model Plot: Log Scale") + ylab("log(Reported Cases + 1)")
ggsave(mle_plot_joint_log, file = "mle_plot_joint_log.png", width = 8, height = 6)
inset_dat <- simdat_test %>%
  dplyr::filter(week > 115)
inset <- plot_fit_full_wribbon(inset_dat) + theme(legend.position = "none", 
                                                  axis.title.x=element_blank(), 
                                                  axis.title.y=element_blank(), 
                                                  axis.text.x = element_text(size = 6),
                                                  axis.text.y = element_text(size = 6))
joint_inset <-
  ggdraw() +
  draw_plot(mle_plot_joint) +
  draw_plot(inset, x = .4, y = .5, width = .55, height = .4)
ggsave(joint_inset, file = "joint_inset.png", width = 8, height = 6)


## adjusted, joint likelihood
coef(joint_mod) <- best_pars
replicate(20, joint_mod %>% pfilter(Np = 5000)) -> reps
epi <- c()
end <- c()
full <- c()
for (i in seq(1,20)) {
  rep <- reps[[i]]
  cond_liks <- rep@cond.logLik
  epi_logliks <- cond_liks[1:232]
  end_logliks <- cond_liks[233:430]
  epi <- c(epi, sum(epi_logliks))
  end <- c(end, sum(end_logliks))
  full <- c(full, logLik(rep))
}
epi_mean <- mean(epi)
end_mean <- mean(end)
full_mean <- logmeanexp(full)

################################################################################
## Forecasting - Figures A2, A3                                               ##
################################################################################
orig_cases <- rep(NA, 1004)
orig_cases[1:430] <- haiti_dat$cases
nsims <- 1000
## no vac
temp1 <- haiti1_joint(vacscen = "id0")
## change times and time zero
fc_mod <- temp1
time(fc_mod) <- seq(1, 1004)
fc_mod@t0 <- 0
## simulate
novac_simdat <- simulate(fc_mod, nsim = nsims, format = "data.frame") %>%
  group_by(`.id`)
fc_simdat <- sim_data(novac_simdat) %>%
  dplyr::rename(est_lo = cases_lo) %>%
  dplyr::rename(est_hi = cases_hi) %>%
  dplyr::rename(est_med = cases_med)
fc_simdat <- fc_simdat %>%
  mutate(period = ifelse(week < 431, "fit", "forecast"))
fc_simdat$orig_cases <- orig_cases
fc_simdat <- fc_simdat %>%
  select(week, est_med, est_lo, est_hi,
         incid_med, incid_lo, incid_hi, orig_cases)
fc_simdat$date <- lubridate::ymd("2010-10-16") + lubridate::weeks(fc_simdat$week)
fc_simdat$vac <- "No Vaccination"
fc_novac <- fc_simdat
forecast <- fc_simdat

## 2 vac
temp2 <- haiti1_joint(vacscen = "id2")
## change times and time zero
fc_mod <- temp2
time(fc_mod) <- seq(1, 1004)
fc_mod@t0 <- 0
depts <- 2
## simulate
vac2_simdat <- simulate(fc_mod, nsim = nsims, format = "data.frame") %>%
  group_by(`.id`)
fc_simdat <- sim_data(vac2_simdat) %>%
  dplyr::rename(est_lo = cases_lo) %>%
  dplyr::rename(est_hi = cases_hi) %>%
  dplyr::rename(est_med = cases_med)
fc_simdat$orig_cases <- NA
fc_simdat$week_end <- NA
fc_simdat <- fc_simdat %>%
  mutate(period = ifelse(week < 431, "endemic", "forecast"))
fc_simdat$orig_cases <- orig_cases
fc_simdat <- fc_simdat %>%
  select(week, est_med, est_lo, est_hi,
         incid_med, incid_lo, incid_hi, orig_cases)
fc_simdat$date <- lubridate::ymd("2010-10-16") + lubridate::weeks(fc_simdat$week)
fc_simdat$vac <- "2 Dept. Vaccination"
fc_2dept <- fc_simdat
forecast <- rbind(forecast, fc_simdat)

## 3 dept
temp3 <- haiti1_joint(vacscen = "id4")
## change times and time zero
fc_mod <- temp3
time(fc_mod) <- seq(1, 1004)
fc_mod@t0 <- 0
depts <- 3
## simulate
vac3_simdat <- simulate(fc_mod, nsim = nsims, format = "data.frame") %>%
  group_by(`.id`)
fc_simdat <- sim_data(vac3_simdat) %>%
  dplyr::rename(est_lo = cases_lo) %>%
  dplyr::rename(est_hi = cases_hi) %>%
  dplyr::rename(est_med = cases_med)
fc_simdat <- fc_simdat %>%
  mutate(period = ifelse(week < 431, "endemic", "forecast"))
fc_simdat$orig_cases <- orig_cases
fc_simdat <- fc_simdat %>%
  select(week, est_med, est_lo, est_hi,
         incid_med, incid_lo, incid_hi, orig_cases)
fc_simdat$date <- lubridate::ymd("2010-10-16") + lubridate::weeks(fc_simdat$week)
fc_simdat$vac <- "3 Dept. Vaccination"
fc_3dept <- fc_simdat
forecast <- rbind(forecast, fc_simdat)

## slow nat
temp4 <- haiti1_joint(vacscen = "id3")
## change times and time zero
fc_mod <- temp4
time(fc_mod) <- seq(1, 1004)
fc_mod@t0 <- 0
depts <- 10
## simulate
slowvac_simdat <- simulate(fc_mod, nsim = nsims, format = "data.frame") %>%
  group_by(`.id`)
fc_simdat <- sim_data(slowvac_simdat) %>%
  dplyr::rename(est_lo = cases_lo) %>%
  dplyr::rename(est_hi = cases_hi) %>%
  dplyr::rename(est_med = cases_med)
fc_simdat <- fc_simdat %>%
  mutate(period = ifelse(week < 431, "endemic", "forecast"))
fc_simdat$orig_cases <- orig_cases
fc_simdat <- fc_simdat %>%
  select(week, est_med, est_lo, est_hi,
         incid_med, incid_lo, incid_hi, orig_cases)
fc_simdat$date <- lubridate::ymd("2010-10-16") + lubridate::weeks(fc_simdat$week)
fc_simdat$vac <- "Slow National Vaccination"
fc_slowvac <- fc_simdat
forecast <- rbind(forecast, fc_simdat)

## fast nat
temp5 <- haiti1_joint(vacscen = "id1")
## change times and time zero
fc_mod <- temp5
time(fc_mod) <- seq(1, 1004)
fc_mod@t0 <- 0
depts <- 10
## simulate
fastvac_simdat <- simulate(fc_mod, nsim = nsims, format = "data.frame") %>%
  group_by(`.id`)
fc_simdat <- sim_data(fastvac_simdat) %>%
  dplyr::rename(est_lo = cases_lo) %>%
  dplyr::rename(est_hi = cases_hi) %>%
  dplyr::rename(est_med = cases_med)
fc_simdat <- fc_simdat %>%
  mutate(period = ifelse(week < 431, "endemic", "forecast"))
fc_simdat$orig_cases <- orig_cases
fc_simdat <- fc_simdat %>%
  select(week, est_med, est_lo, est_hi,
         incid_med, incid_lo, incid_hi, orig_cases)
fc_simdat$date <- lubridate::ymd("2010-10-16") + lubridate::weeks(fc_simdat$week)
fc_simdat$vac <- "Fast National Vaccination"
fc_fastvac <- fc_simdat
forecast <- rbind(forecast, fc_simdat)

## fast high cov
temp6 <- haiti1_joint(vacscen = "id5")
## change times and time zero
fc_mod <- temp6
time(fc_mod) <- seq(1, 1004)
fc_mod@t0 <- 0
depts <- 10
## simulate
hcvac_simdat <- simulate(fc_mod, nsim = nsims, format = "data.frame") %>%
  group_by(`.id`)
fc_simdat <- sim_data(hcvac_simdat) %>%
  dplyr::rename(est_lo = cases_lo) %>%
  dplyr::rename(est_hi = cases_hi) %>%
  dplyr::rename(est_med = cases_med)
fc_simdat <- fc_simdat %>%
  mutate(period = ifelse(week < 431, "endemic", "forecast"))
fc_simdat$orig_cases <- orig_cases
fc_simdat <- fc_simdat %>%
  select(week, est_med, est_lo, est_hi,
         incid_med, incid_lo, incid_hi, orig_cases)
fc_simdat$date <- lubridate::ymd("2010-10-16") + lubridate::weeks(fc_simdat$week)
fc_simdat$vac <- "Fast National, High Coverage Vaccination"
fc_hcvac <- fc_simdat
forecast <- rbind(forecast, fc_simdat)

forecast <- forecast %>%
  mutate(vac = factor(vac, levels = c("No Vaccination",
                                      "2 Dept. Vaccination",
                                      "3 Dept. Vaccination",
                                      "Slow National Vaccination",
                                      "Fast National Vaccination",
                                      "Fast National, High Coverage Vaccination"))) %>%
  filter(week > 430) %>%
  filter(week < 638)

fc_plot <- ggplot(data = forecast, aes(x = date)) +
  geom_line(aes(y = est_med, color = vac), size = 1) +
  geom_ribbon(aes(ymin = est_lo, ymax = est_hi, fill = vac), alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  labs(title = "Forecasts By Vaccination Campaign",
       y = "Reported Cases", x = "Time") +
  facet_wrap(~vac, ncol = 1)

incid_plot <- ggplot(data = forecast, aes(x = date)) +
  geom_line(aes(y = incid_med, color = vac), size = 1) +
  geom_ribbon(aes(ymin = incid_lo, ymax = incid_hi, fill = vac), alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  labs(title = "Forecasts By Vaccination Campaign",
       y = "True Incidence", x = "Time") +
  facet_wrap(~vac, ncol = 1)