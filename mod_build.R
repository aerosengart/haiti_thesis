#==============================================================================#
# Model 1: Model Building                                                      #
#          Code modified from Elizabeth C. Lee, Andrew S. Azman, and Justin    #
#          Lessler at Johns Hopkins Bloomberg School of Public Health          #
#==============================================================================#

####POMP MODEL BUILDING#########################################################
mod_build <- function(dat, my.times = "week", covar.times = "time",
                      my.t0 = 0, covar) {
##----rinit-------------------------------------------------------------------##
  rinit <- Csnippet("
    double frac = pop_0 / (S_0 + E_0 + I_0 + A_0 + R_0);
    S = nearbyint(frac * S_0);
    E = nearbyint(frac * E_0);
    I = nearbyint(frac * I_0);
    A = nearbyint(frac * A_0);
    R = nearbyint(frac * R_0);
    incid = 0;
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
    //cases = rnbinom_mu(tau, rho * incid + 1e-20);
    cases = rnbinom_mu(tau, rho * incid);
    if (cases > 0.0) {
      cases = nearbyint(cases);
    } else {
      cases = 0.0;
    }
  ')

##----dmeasure----------------------------------------------------------------##
  dmeas <- Csnippet('
    if (ISNA(cases)) {
      lik = (give_log) ? 0 : 1;
    } else {
      lik = dnbinom_mu(cases, tau, rho*incid, give_log);
    }
  ')

##----state names-------------------------------------------------------------##
  state_names <- c("S", "E", "I", "A", "R", "incid")

##----parameter names---------------------------------------------------------##
  param_names <- c("rho", "tau", "beta1", "beta2", "beta3",
                   "beta4", "beta5", "beta6", "gamma", "sigma",
                   "theta0", "alpha", "mu", "delta", "nu",
                   "pop_0", "S_0", "E_0", "I_0", "A_0", "R_0")

##----parameter transformations-----------------------------------------------##
  param_trans <- parameter_trans(
    log = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6",
            "tau", "sigma", "gamma", "mu", "delta", "alpha"),
    logit = c("rho", "nu", "theta0"),
    barycentric = c("S_0", "E_0", "I_0", "A_0", "R_0")
  )

##----make pomp model---------------------------------------------------------##
  my_mod <- pomp(
    data = dat,
    times = my.times,
    t0 = my.t0,
    dmeasure = dmeas,
    rmeasure = rmeas,
    rprocess = euler(step.fun = rproc, delta.t = 1/7),
    covar = covariate_table(covar, times = covar.times),
    statenames = state_names,
    paramnames = param_names,
    partrans = param_trans,
    accumvars = c("incid"),
    rinit = rinit
  )
  return(my_mod)
}
