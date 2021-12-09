#==============================================================================#
# Model 1: Helper functions for visualizing simulations                        #
#          Code modified from Elizabeth C. Lee, Andrew S. Azman, and Justin    #
#          Lessler at Johns Hopkins Bloomberg School of Public Health          #
#==============================================================================#

##############################################################################--
####  VISUALIZATION FUNCTIONS FOR ONE FITTING PERIOD                        ####
##############################################################################--
#### plot data for one simulation on a given fitting period
plot_sim_dat <- function(simdat, prof_if, period) {
  if (period == 0) { ## epidemic
    plt <- ggplot(simdat, aes(x = week)) +
      geom_line(aes(y = cases), colour = "blue", alpha = 0.35) +
      geom_point(aes(y = orig), size = 1) +
      annotate("text", x = 100, y = 20000, 
               label = paste0("loglik: ", unlist(prof_if["loglik"]))) + theme_bw() +
      xlim(0, 232) +
      scale_y_continuous("cases", limits = c(0, 30000), oob = scales::squish)
  } else { ## endemic
    plt <- ggplot(simdat, aes(x = week)) +
      geom_line(aes(y = cases), colour = "blue", alpha = 0.35) +
      geom_point(aes(y = orig), size = 1) +
      annotate("text", x = 400, y = 20000, 
               label = paste0("loglik: ", unlist(prof_if["loglik"]))) + theme_bw() +
      xlim(233, 430) +
      scale_y_continuous("cases", limits = c(0, 30000), oob = scales::squish)
  }
  return(plt)
}

#### plot observations with simulated data
##      line - median parameter estimates
##      shaded area - show range of estimates as 2.5 to 97.5 percentiles
plot_fit_wribbon <- function(fitdata) {
  plt <- ggplot(fitdata, aes(x = date)) +
    geom_ribbon(aes(ymin = est_lo, ymax = est_hi), 
                fill = "#458b74", alpha = 0.5) +
    geom_line(aes(y = est_med), size = 1, colour = "#458b74") +
    geom_point(aes(y = orig_cases), colour = "black") +
    scale_x_continuous(breaks = pretty(fitdata$date, n = 20)) +
    theme_bw()
  return(plt)
}

#### plot observations with simulated data
##      line - median number of estimated cases across simulations
plot_fit_noribbon <- function(fitdata) {
  plt <- ggplot(fitdata, aes(x = date)) +
    geom_line(aes(y = est_med), size = 1, colour = "#458b74") +
    geom_point(aes(y = orig_cases), colour = "black") +
    scale_x_continuous(breaks = pretty(fitdata$date, n = 20)) +
    theme_bw()
  return(plt)
}

##############################################################################--
####  VISUALIZATION FUNCTIONS FOR EPIDEMIC AND ENDEMIC PERIODS              ####
##############################################################################--
#### plot observations with simulated data
##      line - median parameter estimates
##      shaded area - show range of estimates as 2.5 to 97.5 percentiles
plot_fit_full_wribbon <- function(fitdata){
  colvec <- c("epidemic" = "#008d23", "endemic" = "#8d006a")
  plt <- ggplot(fitdata, aes(x = date)) +
    ylab("Reported Cases") + 
    geom_ribbon(aes(ymin = est_lo, ymax = est_hi, fill = period), alpha = 0.5) +
    geom_line(aes(y = est_med, colour = period), size = 1) +
    geom_point(aes(y = orig_cases), colour = "black") +
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_continuous(breaks = pretty(fitdata$date, n = 10)) +
    scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))
  return(plt)
}

#### plot observations with simulated data
##      line - median number of estimated cases across simulations
plot_fit_full_noribbon <- function(fitdata){
  colvec <- c("epidemic" = "#008d23", "endemic" = "#8d006a")
  plt <- ggplot(fitdata, aes(x = date)) +
    ylab("Reported Cases") + 
    geom_line(aes(y = est_med, colour = period), size = 1) +
    geom_point(aes(y = orig_cases), colour = "black") +
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_continuous(breaks = pretty(fitdata$date, n = 10)) +
    scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))
  return(plt)
}

##############################################################################--
####  VISUALIZATION FUNCTIONS FOR EPIDEMIC AND ENDEMIC PERIODS              ####
##############################################################################--

## observed fit
plot_fc_case <- function(fdata){
  colvec <- c("epidemic" = "#008d23", "endemic" = "#8d006a", "forecast" = "#8d2300")
  plt <- ggplot(fdata, aes(x = weekdate)) +
    geom_ribbon(aes(ymin = cases_lo, ymax = cases_hi, fill = period), alpha = 0.5) +
    geom_line(aes(y = cases_med, colour = period), size = 1) +
    geom_point(aes(y = orig_cases), colour = "black", alpha = 0.5) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size=11), axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_date(date_breaks = "1 year") +
    labs(x = "Week", y = "Reported, Symptomatic Incidence") +
    scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))
  return(plt)
}

# true forecast
plot_fc_incid <- function(fcdata){
  colvec <- c("epidemic" = "#008d23", "endemic" = "#8d006a", "forecast" = "#8d2300")
  plt <- ggplot(fcdata, aes(x = weekdate)) +
    geom_ribbon(aes(ymin = incid_lo, ymax = incid_hi, fill = period), alpha = 0.5) +
    geom_line(aes(y = incid_med, colour = period), size = 1) +
    geom_point(aes(y = orig_cases), colour = "black", alpha = 0.5) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size=11), axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_date(date_breaks = "1 year") +
    labs(x = "Week", y = "True Cases") +
    scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))
  return(plt)
}


plot_cumvac <- function(vacc_scens){
  
  lastvactime <- vacc_scens %>% dplyr::mutate(any_vacc = s1 + s2 + s3 + s4) %>%
    dplyr::filter(any_vacc > 0) %>%
    dplyr::filter(week == max(week)) %>% dplyr::select(week) %>% unlist %>% unname
  
  deploylabs <- data.frame(scenario = c("s2", "s4", "s3", "s1"), pltlab = c("2 dept", "3 dept", "slow national", "fast national"))
  
  pltdat <- tbl_df(vacc_scens) %>%
    dplyr::filter(week <= lastvactime+1) %>% 
    dplyr::mutate_at(vars(c("s1", "s2", "s3", "s4")), cumsum) %>%
    dplyr::select(weekdate, week, contains("s")) %>%
    gather(scenario, num_vacc, s1:s4) %>%
    dplyr::mutate(scenario = factor(scenario, levels = deploylabs$scenario, labels = deploylabs$pltlab))
  
  
  plt <- ggplot(pltdat, aes(x = weekdate, y = num_vacc, group = scenario)) +
    geom_line(aes(colour = scenario), size = 1) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_colour_viridis_d("deployment strategy", option = "C") +
    scale_y_continuous("Vaccinated People") +
    scale_x_date(date_breaks = "6 months")
  
  return(plt)
  
}

plot_cumvac_frac <- function(vacc_scens, plotdf_fc){
  
  lastvactime <- vacc_scens %>% dplyr::mutate(any_vacc = s1 + s2 + s3 + s4) %>%
    dplyr::filter(any_vacc > 0) %>%
    dplyr::filter(week == max(week)) %>% dplyr::select(week) %>% unlist %>% unname
  
  deploylabs <- data.frame(scenario = c("s2", "s4", "s3", "s1"), pltlab = c("2 dept", "3 dept", "slow national", "fast national"))
  
  pop_df <- plotdf_fc %>% dplyr::select(week, pop_med)
  
  pltdat <- tbl_df(vacc_scens) %>%
    dplyr::filter(week <= lastvactime+1) %>% 
    dplyr::mutate_at(vars(c("s1", "s2", "s3", "s4")), cumsum) %>%
    dplyr::select(weekdate, week, contains("s")) %>%
    gather(scenario, num_vacc, s1:s4) %>%
    left_join(pop_df, by = c("week")) %>%
    dplyr::mutate(scenario = factor(scenario, levels = deploylabs$scenario, labels = deploylabs$pltlab),
                  frac_vacc = num_vacc/pop_med)
  
  
  plt <- ggplot(pltdat, aes(x = weekdate, y = frac_vacc, group = scenario)) +
    geom_line(aes(colour = scenario), size = 1) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_colour_viridis_d("deployment strategy", option = "C") +
    scale_y_continuous("Fraction of Population Vaccinated", limits = c(0,1)) +
    scale_x_date(date_breaks = "6 months")
  
  return(plt)
  
}
