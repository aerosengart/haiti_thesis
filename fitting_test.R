#==============================================================================#
# Model 1: Epidemic Iterated Filtering                                         #
#          Parameter estimation for epidemic period (Oct. 2010 - Mar. 2015)    #
#          Code modified from Elizabeth C. Lee, Andrew S. Azman, and Justin    #
#          Lessler at Johns Hopkins Bloomberg School of Public Health          #
#==============================================================================#
source("utils_clean.R")
source("utils_eval.R")
reload_source()
registerDoParallel(cores=8)

##----settings --- change if desired--------------------------------------------
# filtering and simulation settings
num_iters <- 100 # number of iterations for IF2 (originally 100)
num_parts <- 5000 # number of particles for IF2 (originally 100)
num_sims <- 25 # number of simulations for median calculations (originally 25)

plot_fit_wribbon <- function(fitdata){
  plt <- ggplot(fitdata, aes(x = week)) +
    geom_ribbon(aes(ymin = est_lo, ymax = est_hi),
                fill = "#458b74", alpha = 0.5) +
    geom_line(aes(y = est_med), size = 1, colour = "#458b74") +
    geom_point(aes(y = orig_cases), colour = "black") +
    theme_bw()
  return(plt)
}

plot_fit_noribbon <- function(fitdata){
  plt <- ggplot(fitdata, aes(x = week)) +
    geom_line(aes(y = est_med), size = 1, colour = "#458b74") +
    geom_point(aes(y = orig_cases), colour = "black") +
    theme_bw()
  return(plt)
}

plot_data_epi <- function(df) {
  t <- df %>%
    ungroup %>%
    group_by(week) %>%
    dplyr::summarise(est_med = median(cases),
                     est_lo = quantile(cases, probs = c(.025)),
                     est_hi = quantile(cases, probs = c(.975))) %>%
    dplyr::full_join(haiti.dat.epi %>%
                       dplyr::rename(orig_cases = cases),
                     by = c("week")) %>%
    dplyr::arrange(week)
  return(t)
}

plot_data_end <- function(df) {
  t <- df %>%
    ungroup %>%
    group_by(week_end) %>%
    dplyr::summarise(est_med = median(cases),
                     est_lo = quantile(cases, probs = c(.025)),
                     est_hi = quantile(cases, probs = c(.975))) %>%
    dplyr::full_join(haiti.dat.end %>%
                       dplyr::rename(orig_cases = cases),
                     by = c("week_end")) %>%
    mutate(week = week_end + 232) %>%
    select(-week_end) %>%
    dplyr::arrange(week)
  return(t)
}

####POMP MODEL BUILDING#########################################################
##----rinit-------------------------------------------------------------------##
rinit <- Csnippet("
  double pop = 10911819;
  double frac = pop / (S_0 + E_0 + I_0 + A_0 + R_0);
  S = nearbyint(frac * S_0);
  E = nearbyint(frac * E_0);
  I = nearbyint(frac * I_0);
  A = nearbyint(frac * A_0);
  R = nearbyint(frac * R_0);
  incid = I; // initialize as entrants into I state
")

##----rprocess----------------------------------------------------------------##
rproc <- Csnippet("
  // transition rates
  double S_rate[2], E_rate[3], I_rate[2], A_rate[2], R_rate[2];
  // transition numbers
  double S_trans[2], E_trans[3], I_trans[2], A_trans[2], R_trans[2];
  
  // seasonality
  double beta = beta1*seas1 + beta2*seas2 + 
                beta3*seas3 + beta4*seas4 + 
                beta5*seas5 + beta6*seas6;
  
  // population and births
  int pop = S + E + I + A + R;
  int births = rpois(mu * pop * dt);
  
  // force of infection
  double foi = pow(I, nu) * (beta / pop);
  
  // transition rate calculations
  S_rate[0] = foi;  // S -> E
  S_rate[1] = delta; // S -> death

  E_rate[0] = sigma * (1 - theta0); // E -> I
  E_rate[1] = sigma * theta0; // E -> A
  E_rate[2] = delta; // E -> death

  I_rate[0] = gamma; // I -> R
  I_rate[1] = delta; // I -> death

  A_rate[0] = gamma; // A -> R
  A_rate[1] = delta; // A -> death

  R_rate[0] = alpha; // R -> S 
  R_rate[1] = delta; // R -> death 
  
  // transition numbers
  reulermultinom(2, S, &S_rate[0], dt, &S_trans[0]);
  reulermultinom(3, E, &E_rate[0], dt, &E_trans[0]);
  reulermultinom(2, I, &I_rate[0], dt, &I_trans[0]);
  reulermultinom(2, A, &A_rate[0], dt, &A_trans[0]);
  reulermultinom(2, R, &R_rate[0], dt, &R_trans[0]);
  
  // states
  // S += - (S to E) - (S to death) + (R to S) + births
  S += -S_trans[0] - S_trans[1] + R_trans[0] + births;
  // E += - (E to I) - (E to A) - (E to death) + (S to E)
  E += -E_trans[0] - E_trans[1] - E_trans[2] + S_trans[0];
  // I += - (I to R) - (I to death) + (E to I)
  I += -I_trans[0] - I_trans[1] + E_trans[0];
  // A += - (A to R) - (A to death) + (E to R)
  A += -A_trans[0] - A_trans[1] + E_trans[1];
  // R += - (R to S) - (R to death) + (I to R) + (A to R)
  R += -R_trans[0] - R_trans[1] + I_trans[0] + A_trans[0];
  
  // incidence is cumulative entries into I state
  incid += E_trans[0]; 
")

##----rmeasure----------------------------------------------------------------##
rmeas <- Csnippet("
  cases = rnbinom_mu(tau, rho*incid);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } 
  else {
    cases = 0.0;
  }
")

##----dmeasure----------------------------------------------------------------##
dmeas <- Csnippet("
  lik = dnbinom_mu(cases, tau, rho*incid, give_log);
")

##----skeleton----------------------------------------------------------------##
sim.skel <- Csnippet('
  // transition rates
  double S_rate[2]; 
  double E_rate[3];
  double I_rate[2];
  double A_rate[2];
  double R_rate[2];

  // transition terms
  double S_term[2]; 
  double E_term[3];
  double I_term[2];
  double A_term[2];
  double R_term[2];

  // some population demonitors
  int pop = S + E + I + A + R;
  double births = mu * pop;

  // make seasonal beta term for current time
  double beta = beta1*seas1 + beta2*seas2 + 
                beta3*seas3 + beta4*seas4 + 
                beta5*seas5 + beta6*seas6;
  double foi = pow(I, nu) * beta / pop;

  //compute the rates for all the transitions
  S_rate[0]= foi;  // S -> E
  S_rate[1]= delta; // S -> dead

  E_rate[0]= sigma*(1-theta0); // E -> I
  E_rate[1]= sigma*theta0; // E -> A
  E_rate[2]= delta; // E -> dead

  I_rate[0]= gamma; // I -> R
  I_rate[1]= delta;

  A_rate[0]= gamma; // A -> R
  A_rate[1]= delta; // A -> dead

  R_rate[0]= alpha; // R                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         -> S exp. waning immunity from natural infection
  R_rate[1]= delta; // R -> dead

  // compute the transition terms
  for(int i=0; i < 2; i++) {
    double term = S_rate[i] * S;
    S_term[i] = term;
  }
  for(int i=0; i < 3; i++) {
    double term = E_rate[i] * E;
    E_term[i] = term;
  }
  for(int i=0; i < 2; i++) {
    double term = I_rate[i] * I;
    I_term[i] = term;
  }
  for(int i=0; i < 2; i++) {
    double term = A_rate[i] * A;
    A_term[i] = term;
  }
  for(int i=0; i < 2; i++) {
    double term = R_rate[i] * R;
    R_term[i] = term;
  }

  // balance the equations
  DS = -S_term[0] - S_term[1] + R_term[0] + births; 
  DE = -E_term[0] - E_term[1] - E_term[2] + S_term[0];
  DI = -I_term[0] - I_term[1] + E_term[0];
  DA = -A_term[0] - A_term[1] + E_term[1];
  DR = -R_term[0] - R_term[1] + I_term[0] + A_term[0];
  Dincid = E_term[0]; // incidence is cumulative entries into I state
')

##----state names-------------------------------------------------------------##
state_names <- c("S", "E", "I", "A", "R", "incid")

##----parameter names---------------------------------------------------------##
param_names <- c("rho", "tau", "beta1", "beta2", "beta3",
                 "beta4", "beta5", "beta6", "gamma", "sigma",
                 "theta0", "alpha", "mu", "delta", "nu",
                 "S_0", "E_0", "I_0", "A_0", "R_0")

##----parameter transformations-----------------------------------------------##
param_trans <- parameter_trans(
  log = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6",
          "tau", "sigma", "gamma", "mu", "delta", "alpha"),
  logit = c("rho", "nu", "theta0"),
  log_barycentric = c("S_0", "E_0", "I_0", "A_0", "R_0")
)

##----make pomp model---------------------------------------------------------##
build.epi.mod <- function(pop, ##  = get.haiti.pop()
                          dat, ##  = get.mspp.agg.data()
                          my.times = "week", 
                          my.t0 = 0,
                          covar,
                          pars){
  my.mod <- pomp(
    data = dat,
    times = my.times,
    t0 = my.t0,
    dmeasure = dmeas,
    rmeasure = rmeas,
    rprocess = euler(step.fun=rproc, delta.t=1/7),
    skeleton = vectorfield(sim.skel),
    covar = covariate_table(covar, times = "time"),
    partrans = param_trans,
    statenames = state_names,
    params = pars,
    paramnames = param_names,
    accumvars = c("incid"),
    rinit = rinit
  )
  return(my.mod)
}

########
pop.haiti <- get.haiti.pop()
haiti.dat <- get.mspp.agg.data() 
covartab <- make.covartab(0, nrow(haiti.dat.epi) + 1, byt = 1, 
                          degree = 6, nbasis = 6, per = 52.14)
haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= 232)
# parameters from supplement
supp_e0 <- 405/pop.haiti; supp_i0 <- 518/pop.haiti
supp_params <- c(rho = 0.34, tau = 4.0, beta1 = 5.8, beta2 = 3.5, beta3 = 4.4,
            beta4 = 2.8, beta5 = 6.2, beta6 = 2.4, nu = 0.96, gamma = 7/2, 
            sigma = 7/1.4, theta0 = 0, alpha = 7/920, 
            mu = ((1+22.6/1000)^(1/52.14))-1, delta = ((1+7.5/1000)^(1/52.14))-1,
            S_0 = 1 - supp_e0 - supp_i0, E_0 = supp_e0, I_0 = supp_i0, A_0 = 0.0, R_0 = 0.0)
# parameters from MIF2 AER
aer_params <- c(rho = 0.7476200871387470, 
                 tau = 2.8147689942898600,
                 beta1 = 2.656989965595180,
                 beta2 = 6.52646475897578,
                 beta3 = 1.001065528981810,
                 beta4 = 6.212120297172650, 
                 beta5 = 2.0735618284143900, 
                 beta6 = 4.982623401556510,
                 nu = 0.951038893143546, 
                 gamma = 7/2, sigma = 7/1.4, 
                 theta0 = 0, alpha = 7/920, 
                 mu = ((1+22.6/1000)^(1/52.14))-1, delta = ((1+7.5/1000)^(1/52.14))-1,
                 S_0 = 0.9647493813555620, 
                E_0 = 3.60213029985934E-05, 
                I_0 = 4.57985138124973E-05,
                A_0 = 0.0, 
                R_0 = 0.03516879882762670)
# parameters from MIF2 AER
# build models
supp_mod_epi <- build.epi.mod(pop = pop.haiti, dat = haiti.dat.epi,
                               covar = covartab, pars = supp_params)
aer_mod_epi <- build.epi.mod(pop = pop.haiti, dat = haiti.dat.epi,
                              covar = covartab, pars = aer_params)

# simulation
supp_sim_epi <- simulate(supp_mod_epi, params = supp_params, nsim = 10, format = "data.frame")
aer_sim_epi <- simulate(aer_mod_epi, params = aer_params, nsim = 10, format = "data.frame")

# simulation visual
plot1 <- supp_sim_epi %>% 
  plot_data_epi() %>%
  plot_fit_wribbon()
ggsave("testing/simulation_epi_supplement_params.pdf", plot1)

plot2 <- aer_sim_epi %>% 
  plot_data_epi() %>%
  plot_fit_wribbon()
ggsave("testing/simulation_epi_aer_params.pdf", plot2)

# filter again
foreach(i = 1:10, .combine = c) %dopar% {
  supp_mod_epi %>%
    pfilter(params = supp_params, Np = 5000)
} -> pf10_supp
foreach(i = 1:10, .combine = c) %dopar% {
  aer_mod_epi %>%
    pfilter(params = aer_params, Np = 5000)
} -> pf10_aer

L_pf10_supp <- pf10_supp %>%
  logLik() %>%
  logmeanexp(se = TRUE) # getting an unbiased likelihood estimate of -2300
L_pf10_aer <- pf10_aer %>%
  logLik() %>%
  logmeanexp(se = TRUE) # getting an unbiased likelihood estimate of -2270

# save likelihood results for this point
pf10_supp[[1]] %>% coef() %>% bind_rows() %>%
  bind_cols(loglik = L_pf10_supp[1], loglik.se = L_pf10_supp[2]) %>%
  write_csv("testing/pf10_supp_fit_test_params.csv")
pf10_aer[[1]] %>% coef() %>% bind_rows() %>%
  bind_cols(loglik = L_pf10_aer[1], loglik.se = L_pf10_aer[2]) %>%
  write_csv("testing/pf10_aer_fit_test_params.csv")

# local search of likelihood surface
foreach(i = 1:10, .combine = c) %dopar% {
  supp_mod_epi %>%
    mif2(
      params = supp_params,
      Np = 1000, Nmif = 50,
      cooling.fraction.50 = 0.5,
      rw.sd = rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02, beta4 = 0.02,
                    beta5 = 0.02, beta6 = 0.02, tau = 0.005, rho = 0.01, nu = 0.04,
                    E_0 = ivp(2E-7), I_0 = ivp(2E-7))
    )
} -> mifs_local_supp
foreach(i = 1:10, .combine = c) %dopar% {
  aer_mod_epi %>%
    mif2(
      params = aer_params,
      Np = 1000, Nmif = 50,
      cooling.fraction.50 = 0.5,
      rw.sd = rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02, beta4 = 0.02,
                    beta5 = 0.02, beta6 = 0.02, tau = 0.005, rho = 0.01, nu = 0.04,
                    E_0 = ivp(2E-7), I_0 = ivp(2E-7))
    )
} -> mifs_local_aer

# diagnostics
pdf("testing/mifs_localsearch_supp.pdf")
plot(mifs_local_supp)
dev.off()
pdf("testing/mifs_localsearch_aer.pdf")
plot(mifs_local_aer)
dev.off()

# estimating likelihood from mif2
foreach(mf = mifs_local_supp, .combine = rbind) %dopar% {
  evals <- replicate(10, logLik(pfilter(mf, Np = 5000)))
  ll <- logmeanexp(evals, se = TRUE)
  mf %>% coef() %>% bind_rows() %>%
    bind_cols(loglik = ll[1], loglik.se = ll[2])
} -> results_supp

foreach(mf = mifs_local_aer, .combine = rbind) %dopar% {
  evals <- replicate(10, logLik(pfilter(mf, Np = 5000)))
  ll <- logmeanexp(evals, se = TRUE)
  mf %>% coef() %>% bind_rows() %>%
    bind_cols(loglik = ll[1], loglik.se = ll[2])
} -> results_aer

# look at geometry of likelihood surface in neighborhood of point estimate
pairs(~loglik + beta1 + rho + nu + tau, data = results_supp, pch = 16)

# add this to csv
read_csv("testing/pf10_supp_fit_test_params.csv") %>%
  bind_rows(results_supp) %>%
  arrange(-loglik) %>%
  write_csv("testing/supp_fit_test_params.csv")
read_csv("testing/pf10_aer_fit_test_params.csv") %>%
  bind_rows(results_aer) %>%
  arrange(-loglik) %>%
  write_csv("testing/aer_fit_test_params.csv")

# global search for MLE
guesses <- runif_design(
  lower = c(rho = 1E-8, tau = 1, beta1 = 1E-9, beta2 = 1E-9, 
            beta3 = 1E-9, beta4 = 1E-9, beta5 = 1E-9, beta6 = 1E-9,
            nu = 0.95, E_0 = 1/pop.haiti, I_0 = 1/pop.haiti),
  upper = c(rho = 1, tau = 20, beta1 = 10, beta2 = 10, 
            beta3 = 10, beta4 = 10, beta5 = 10, beta6 = 10,
            nu = 1, E_0 = 500/pop.haiti, I_0 = 500/pop.haiti),
  nseq = 300
)

# Fixed parameters
fix_params <- c(pop_0 = pop.haiti, 
                mu = ((1+22.6/1000)^(1/52.14))-1, 
                delta = ((1+7.5/1000)^(1/52.14))-1,
                sigma = 7/1.4, 
                gamma = 7/2,
                theta0 = 0,
                alpha = 7/2920,
                S_0 = unname(coef(aer_mod_epi, "S_0")),
                A_0 = 0.0,
                R_0 = unname(coef(aer_mod_epi, "R_0")))

mf1_supp <- mifs_local_supp[[1]]
mf1_aer <- mifs_local_aer[[1]]
mf1_aer_res <- vector(length = 300)
i <- 0

# do search
#### one search from each of 300 starting values in guesses
#### each search has initial run of 50 IF2 iterations then another 100 iterations
#### use particle filter to evaluate likelihood
#### return endpoint of search with likelihood estimate and standard error
registerDoRNG(787452981)
foreach(guess = iter(guesses, "row"), .combine = rbind) %dopar% {
  i <- i + 1
  mf1_aer %>%
    mif2(params = c(unlist(guess), fix_params), partrans = param_trans) -> mf
  evals <- replicate(10, logLik(pfilter(mf, Np = 1000)))
  ll <- logmeanexp(evals, se = TRUE)
  mf %>% coef() %>% bind_rows() %>%
    bind_cols(loglik = ll[1], loglik.se = ll[2]) -> t
  mf1_aer_res[i] <- t
} -> glob_results_aer

# look at geometry of likelihood surface
read_csv("fit_test_params.csv") %>%
  filter(loglik > max(loglik) - 50) %>%
  bind_rows(guesses) %>%
  mutate(type = if_else(if.na(loglik), "guess", "result")) %>%
  arrange(type) -> all
# starting values in grey, IF2 estimates in red
pairs(~loglik + beta1 + rho + nu + tau, data = all,
      col = ifelse(all$type == "guess", grey(0.5), "red"), pch = 16)

# poor man's profile
all %>%
  filter(type == "result") %>%
  filter(loglik > max(loglik) - 10) %>%
  ggplot(aes(x = rho, y = loglik)) +
  geom_point() +
  labs(
    x = expression("rho"),
    title = "poor man's profile likelihood - rho"
  )

#### ARMA benchmark
a201 <- arima(df$casesSum, order = c(2, 0, 1))
log_a201 <- arima(log(df$casesSum), order = c(2, 0, 1))


##################################################################################################################
####  ENDEMIC FITTING  ####
##################################################################################################################
##----rinit-------------------------------------------------------------------##
rinit <- Csnippet("
  double pop = pop_0;
  double frac = pop / (S_0 + E_0 + I_0 + A_0 + R_0);
  S = nearbyint(frac * S_0);
  E = nearbyint(frac * E_0);
  I = nearbyint(frac * I_0);
  A = nearbyint(frac * A_0);
  R = nearbyint(frac * R_0);
  incid = incid_0;
")

##----rprocess----------------------------------------------------------------##
rproc <- Csnippet('
  // transition rates
  double Srate[2]; 
  double Erate[3];
  double Irate[2];
  double Arate[2];
  double Rrate[2];

  // transition numbers
  double Strans[2]; 
  double Etrans[3];
  double Itrans[2];
  double Atrans[2];
  double Rtrans[2];

  // some population demonitors
  int pop = S + E + I + A + R;
  int births = rpois(mu*pop*dt);

  // make seasonal beta term for current time
  double mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + 
                  beta4*seas4 + beta5*seas5 + beta6*seas6;
  double foi = pow(I, nu) * mybeta / pop;  

  //compute the rates for all the transitions
  Srate[0]= foi;  // S -> E
  Srate[1]= delta; // S -> death

  Erate[0]= sigma*(1-theta0); // E -> I
  Erate[1]= sigma*theta0; // E -> A
  Erate[2]= delta; // E -> death

  Irate[0]= gamma; // I -> R
  Irate[1]= delta; // I -> death

  Arate[0]= gamma; // A -> R
  Arate[1]= delta; // A -> death

  Rrate[0]= alpha; // R -> S
  Rrate[1]= delta; // R -> death

  // compute the transition numbers
  reulermultinom(2, S, &Srate[0], dt, &Strans[0]);
  reulermultinom(3, E, &Erate[0], dt, &Etrans[0]);
  reulermultinom(2, I, &Irate[0], dt, &Itrans[0]);
  reulermultinom(2, A, &Arate[0], dt, &Atrans[0]);
  reulermultinom(2, R, &Rrate[0], dt, &Rtrans[0]);

  // S += - (S to E) - (S to death) + (R to S) + births
  S += -Strans[0] - Strans[1] + Rtrans[0] + births;
  // E += - (E to I) - (E to A) - (E to death) + (S to E)
  E += -Etrans[0] - Etrans[1] - Etrans[2] + Strans[0];
  // I += - (I to R) - (I to death) + (E to I)
  I += -Itrans[0] - Itrans[1] + Etrans[0];
  // A += - (A to R) - (A to death) + (E to A)
  A += -Atrans[0] - Atrans[1] + Etrans[1];
  // R += - (R to S) - (R to death) + (I to R) + (A to R)
  R += -Rtrans[0] - Rtrans[1] + Itrans[0] + Atrans[0];

  // incidence is cumulative entries into I state
  incid += Etrans[0];
')

##----rmeasure----------------------------------------------------------------##
rmeas <- Csnippet('
  cases = rnbinom_mu(tau, rho * incid + 1e-20);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
')

##----dmeasure----------------------------------------------------------------##
dmeas <- Csnippet('
  lik = dnbinom_mu(cases, tau, rho * incid + 1e-20, give_log);
')

##----skeleton----------------------------------------------------------------##
sim.skel <- Csnippet('
  // transition rates
  double Srate[2]; 
  double Erate[3];
  double Irate[2];
  double Arate[2];
  double Rrate[2];

  // transition terms
  double Sterm[2]; 
  double Eterm[3];
  double Iterm[2];
  double Aterm[2];
  double Rterm[2];

  // some population demonitors
  int pop = S + E + I + A + R;
  double births = mu * pop;

  // make seasonal beta term for current time
  double mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + 
                  beta4*seas4 + beta5*seas5 + beta6*seas6;
  double foi = pow(I, nu)*mybeta/pop;  

  //compute the rates for all the transitions
  Srate[0]= foi;  // S -> E
  Srate[1]= delta; // S -> death

  Erate[0]= sigma*(1-theta0); // E -> I
  Erate[1]= sigma*theta0; // E -> A
  Erate[2]= delta; // E -> death

  Irate[0]= gamma; // I -> R
  Irate[1]= delta; // I -> death

  Arate[0]= gamma; // A -> R
  Arate[1]= delta; // A -> death

  Rrate[0]= alpha; // R -> S
  Rrate[1]= delta; // R -> death

  // compute the transition terms
  for(int i=0; i < 2; i++) {
    double term = Srate[i] * S;
    Sterm[i] = term;
  }
  for(int i=0; i < 3; i++) {
    double term = Erate[i] * E;
    Eterm[i] = term;
  }
  for(int i=0; i < 2; i++) {
    double term = Irate[i]*I;
    Iterm[i] = term;
  }
  for(int i=0; i < 2; i++) {
    double term = Arate[i] * A;
    Aterm[i] = term;
  }
  for(int i=0; i < 2; i++) {
    double term = Rrate[i] * R;
    Rterm[i] = term;
  }

  // balance the equations
  DS = -Sterm[0] - Sterm[1] + Rterm[0] + births; 
  DE = -Eterm[0] - Eterm[1] - Eterm[2] + Sterm[0];
  DI = -Iterm[0] - Iterm[1] + Eterm[0];
  DA = -Aterm[0] - Aterm[1] + Eterm[1];
  DR = -Rterm[0] - Rterm[1] + Iterm[0] + Aterm[0];
  Dincid = Eterm[0]; // incidence is cumulative entries into I state
')

##----state names-------------------------------------------------------------##
state_names <- c("S", "E", "I", "A", "R", "incid")

##----parameter names---------------------------------------------------------##
param_names <- c("rho", "tau", "beta1", "beta2", "beta3",
                 "beta4", "beta5", "beta6", "gamma", "sigma",
                 "theta0", "alpha", "mu", "delta", "nu", "incid_0",
                 "pop_0", "S_0", "E_0", "I_0", "A_0", "R_0")

##----parameter transformations-----------------------------------------------##
param_trans <- parameter_trans(
  log = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6",
          "tau", "sigma", "gamma", "mu", "delta", "alpha"),
  logit = c("rho", "nu", "theta0"),
  log_barycentric = c("S_0", "E_0", "I_0", "A_0", "R_0") ## need?
)

##----make pomp model---------------------------------------------------------##
build.end.mod <- function(dat, ##  = get.mspp.agg.data()
                          my.times = "week_end", 
                          covar.times = "time_end", my.t0 = 0, covar, pars) {
  # Get initial parameter values
  #pars <- pars[1, ]
  
  my.mod <- pomp(
    data = dat,
    times = my.times,
    t0 = my.t0,
    dmeasure = dmeas,
    rmeasure = rmeas,
    rprocess = euler(step.fun=rproc, delta.t=1/7),
    skeleton = vectorfield(sim.skel),
    covar = covariate_table(covar, times = covar.times),
    statenames = state_names,
    paramnames = param_names,
    partrans = param_trans,
    params = pars,
    accumvars = c("incid"),
    rinit = rinit
  )
  return(my.mod)
}

########
haiti.dat.end <- haiti.dat %>% dplyr::filter(week > 232) %>%
  dplyr::mutate(week_end = seq_along(week)) %>%
  select(-week)
covartab.all <- make.covartab(0, nrow(haiti.dat) + 1, byt = 1, degree = 6, 
                              nbasis = 6, per = 52.14)
covartab <- covartab.all[(episettings$nweeks + 1):nrow(covartab.all),] %>%
  dplyr::mutate(time_end = seq_along(time)-1) %>%
  dplyr::select(time_end, contains("seas"))
num_betas <- ncol(covartab.all) - 1

# add/update parameters
supp_params["pop_0"] <- fit_states_supp[232, "pop_med"]
supp_params["incid_0"] <- fit_states_supp[232, "incid_med"]
supp_params["S_0"] <- fit_states_supp[232, "S_med"] / supp_params["pop_0"]
supp_params["E_0"] <- fit_states_supp[232, "E_med"] / supp_params["pop_0"]
supp_params["I_0"] <- fit_states_supp[232, "I_med"] / supp_params["pop_0"]
supp_params["R_0"] <- fit_states_supp[232, "R_med"] / supp_params["pop_0"]

aer_params["pop_0"] <- fit_states_aer[232, "pop_med"]
aer_params["incid_0"] <- fit_states_aer[232, "incid_med"]
aer_params["S_0"] <- fit_states_aer[232, "S_med"] / aer_params["pop_0"]
aer_params["E_0"] <- fit_states_aer[232, "E_med"] / aer_params["pop_0"]
aer_params["I_0"] <- fit_states_aer[232, "I_med"] / aer_params["pop_0"]
aer_params["R_0"] <- fit_states_aer[232, "R_med"] / aer_params["pop_0"]

# build models
supp_mod_end <- build.end.mod(dat = haiti.dat.end,
                              covar = covartab, pars = supp_params)
aer_mod_end <- build.end.mod(dat = haiti.dat.end,
                             covar = covartab, pars = aer_params)

# initial simulation just using epidemic parameter estimates
supp_sim_end <- simulate(supp_mod_end, params = supp_params, nsim = 10, format = "data.frame") 
aer_sim_end <- simulate(aer_mod_end, params = aer_params, nsim = 10, format = "data.frame")

# visualize
plot1 <- supp_sim_end %>% 
  plot_data_end() %>%
  plot_fit_wribbon()
ggsave("testing/simulation_end_supplement_epi_params.pdf", plot1)

plot2 <- aer_sim_end %>% 
  plot_data_end() %>%
  plot_fit_wribbon()
ggsave("testing/simulation_epi_aer_epi_params.pdf", plot2)

# filter again
foreach(i = 1:10, .combine = c) %dopar% {
  supp_mod_end %>%
    pfilter(params = supp_params, Np = 5000)
} -> pf10_supp
foreach(i = 1:10, .combine = c) %dopar% {
  aer_mod_end %>%
    pfilter(params = aer_params, Np = 5000)
} -> pf10_aer

L_pf10_supp <- pf10_supp %>%
  logLik() %>%
  logmeanexp(se = TRUE) # getting an unbiased likelihood estimate of -2300
L_pf10_aer <- pf10_aer %>%
  logLik() %>%
  logmeanexp(se = TRUE) # getting an unbiased likelihood estimate of -1500

# save likelihood results for this point
pf10_supp[[1]] %>% coef() %>% bind_rows() %>%
  bind_cols(loglik = L_pf10_supp[1], loglik.se = L_pf10_supp[2]) %>%
  write_csv("testing/pf10_supp_fit_epi_params.csv")
pf10_aer[[1]] %>% coef() %>% bind_rows() %>%
  bind_cols(loglik = L_pf10_aer[1], loglik.se = L_pf10_aer[2]) %>%
  write_csv("testing/pf10_aer_fit_epi_params.csv")

####
# parameters from supplement reports
coef(supp_mod_end, "rho") <- 0.34
coef(supp_mod_end, "tau") <- 4.1
coef(supp_mod_end, "beta1") <- 5.4
coef(supp_mod_end, "beta2") <- 3.3
coef(supp_mod_end, "beta3") <- 4.3
coef(supp_mod_end, "beta4") <- 3.5
coef(supp_mod_end, "beta5") <- 5.0
coef(supp_mod_end, "beta6") <- 2.7
coef(supp_mod_end, "nu") <- 0.96

# parameters from aer MIF2 -> parid = 297
coef(aer_mod_end, "rho") <- 0.7431591577641590
coef(aer_mod_end, "tau") <- 1.1951978386882100
coef(aer_mod_end, "beta1") <- 3.940021059877410
coef(aer_mod_end, "beta2") <- 3.5387601023490800
coef(aer_mod_end, "beta3") <- 3.3888485776216100
coef(aer_mod_end, "beta4") <- 4.003675189051050
coef(aer_mod_end, "beta5") <- 3.6848009538591000
coef(aer_mod_end, "beta6") <- 3.222459132291050
coef(aer_mod_end, "nu") <- 0.9374121210475350

# simulation from calibrations
supp_sim_end <- simulate(supp_mod_end, params = supp_params, nsim = 10, format = "data.frame")
aer_sim_end <- simulate(aer_mod_end, params = aer_params, nsim = 10, format = "data.frame")

# visualize
plot1 <- supp_sim_end %>% 
  plot_data_end() %>%
  plot_fit_wribbon()
ggsave("testing/simulation_end_supplement_params.pdf", plot1)

plot2 <- aer_sim_end %>% 
  plot_data_end() %>%
  plot_fit_wribbon()
ggsave("testing/simulation_epi_aer_params.pdf", plot2)

# filter again
foreach(i = 1:10, .combine = c) %dopar% {
  supp_mod_end %>%
    pfilter(params = supp_params, Np = 5000)
} -> pf10_supp
foreach(i = 1:10, .combine = c) %dopar% {
  aer_mod_end %>%
    pfilter(params = aer_params, Np = 5000)
} -> pf10_aer

L_pf10_supp <- pf10_supp %>%
  logLik() %>%
  logmeanexp(se = TRUE) # getting an unbiased likelihood estimate of -2300
L_pf10_aer <- pf10_aer %>%
  logLik() %>%
  logmeanexp(se = TRUE) # getting an unbiased likelihood estimate of -1560

# save likelihood results for this point
pf10_supp[[1]] %>% coef() %>% bind_rows() %>%
  bind_cols(loglik = L_pf10_supp[1], loglik.se = L_pf10_supp[2]) %>%
  write_csv("testing/pf10_supp_fit_test_params.csv")
pf10_aer[[1]] %>% coef() %>% bind_rows() %>%
  bind_cols(loglik = L_pf10_aer[1], loglik.se = L_pf10_aer[2]) %>%
  write_csv("testing/pf10_aer_fit_test_params.csv")

# local search of likelihood surface
registerDoRNG(767364278)
foreach(i = 1:20, .combine = c) %dopar% {
  supp_mod_end %>%
    mif2(
      params = supp_params,
      Np = 1000, Nmif = 50,
      cooling.fraction.50 = 0.5,
      rw.sd = rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02, beta4 = 0.02,
                    beta5 = 0.02, beta6 = 0.02, tau = 0.005, rho = 0.01, nu = 0.04,
                    E_0 = ivp(2E-7), I_0 = ivp(2E-7))
    )
} -> mifs_local_supp
foreach(i = 1:20, .combine = c) %dopar% {
  aer_mod_end %>%
    mif2(
      params = aer_params,
      Np = 1000, Nmif = 50,
      cooling.fraction.50 = 0.5,
      rw.sd = rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02, beta4 = 0.02,
                    beta5 = 0.02, beta6 = 0.02, tau = 0.005, rho = 0.01, nu = 0.04,
                    E_0 = ivp(2E-7), I_0 = ivp(2E-7))
    )
} -> mifs_local_aer

# diagnostics
pdf("testing/mifs_localsearch_supp.pdf")
plot(mifs_local_supp)
dev.off()
pdf("testing/mifs_localsearch_aer.pdf")
plot(mifs_local_aer)
dev.off()

# estimating likelihood from mif2
foreach(mf = mifs_local_supp, .combine = rbind) %dopar% {
  evals <- replicate(10, logLik(pfilter(mf, Np = 5000)))
  ll <- logmeanexp(evals, se = TRUE)
  mf %>% coef() %>% bind_rows() %>%
    bind_cols(loglik = ll[1], loglik.se = ll[2])
} -> results_supp

foreach(mf = mifs_local_aer, .combine = rbind) %dopar% {
  evals <- replicate(10, logLik(pfilter(mf, Np = 5000)))
  ll <- logmeanexp(evals, se = TRUE)
  mf %>% coef() %>% bind_rows() %>%
    bind_cols(loglik = ll[1], loglik.se = ll[2])
} -> results_aer

# look at geometry of likelihood surface in neighborhood of point estimate
pairs(~loglik + beta1 + rho + nu + tau, data = results_supp, pch = 16)

# add this to csv
read_csv("testing/pf10_supp_fit_test_params.csv") %>%
  bind_rows(results_supp) %>%
  arrange(-loglik) %>%
  write_csv("testing/supp_fit_test_params.csv")
read_csv("testing/pf10_aer_fit_test_params.csv") %>%
  bind_rows(results_aer) %>%
  arrange(-loglik) %>%
  write_csv("testing/aer_fit_test_params.csv")

# global search for MLE
guesses <- runif_design(
  lower = c(rho = 1E-8, tau = 1, beta1 = 1E-9, beta2 = 1E-9, 
            beta3 = 1E-9, beta4 = 1E-9, beta5 = 1E-9, beta6 = 1E-9,
            nu = 0.95, E_0 = 1/pop.haiti, I_0 = 1/pop.haiti),
  upper = c(rho = 1, tau = 20, beta1 = 10, beta2 = 10, 
            beta3 = 10, beta4 = 10, beta5 = 10, beta6 = 10,
            nu = 1, E_0 = 500/pop.haiti, I_0 = 500/pop.haiti),
  nseq = 300
)

mf1_supp <- mifs_local_supp[[1]]
mf1_aer <- mifs_local_aer[[1]]
mf1_aer_res <- vector(length = 300)
i <- 0
