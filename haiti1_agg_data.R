#' Clean haiti cases data for Model 1

haiti1_agg_data <- function(){
  load("data/haiti_case_data.rda")
  allDat <- haiti_case_data

  splitDate <- strsplit(allDat$date_sat_orig, "-")
  data.table::setattr(splitDate[[1]], 'names', c("year", "month", "day"))
  dateDf <- tibble::as_tibble(as.data.frame(do.call(rbind, splitDate))) %>%
    dplyr::mutate(month = as.character(month)) %>%
    dplyr::mutate(day = as.character(day)) %>%
    dplyr::mutate(year = as.character(year)) %>%
    dplyr::mutate(month = ifelse(nchar(month) == 1, paste0("0", month), month)) %>%
    dplyr::mutate(day = ifelse(nchar(day) == 1, paste0("0", day), day)) %>%
    dplyr::mutate(date_sat = as.Date(paste(year, month, day, sep = "-"), origin = "1900-01-01"))

  fullDateVec <- data.frame(date_sat = seq(min(dateDf$date_sat), max(dateDf$date_sat), by=7))

  cleanDat <- allDat %>%
    dplyr::mutate(date_sat = as.Date(dateDf$date_sat, origin = "1900-01-01")) %>%
    dplyr::select(-date_sat_orig) %>%
    dplyr::full_join(fullDateVec, by = c("date_sat")) %>%
    dplyr::arrange(date_sat) %>%
    dplyr::mutate(week = seq_along(date_sat)) %>%
    tidyr::gather(department, cases, Artibonite:Sud_Est)

  aggDat <- cleanDat %>%
    dplyr::group_by(week) %>% ## "day" or "week"
    dplyr::summarise(cases = sum(cases))

  return(aggDat)
}
