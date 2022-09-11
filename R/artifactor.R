# Scripts on this file are dedicated to performing M-CAPA analysis
# on EEG data.

library(anomaly)
library(methods)
library(tidyverse)
source("R/eeg.R")

source("R/cln.R")
source("R/analysis.r")


#' Perform M-CAPA analysis on EEG data in a given range of seconds
#' and return results filtered by anomaly strength.
#'
#' @param eeg An eeg object
#' @param s First second of timespan to analyze
#' @param e Last second of timespan to analyze
#' @param res Resolution at which to perform analysis
#' @param alpha Threshold of strength significance for collective anomalies
#' @param beta Threshold of strength significance for point anomalies
#' @param thresh How many seconds collective anomaly n must be from
#' collective anomaly (n - 1) to consider them part of a same cluster?

#' @return An analysis object
analyze <- function(eeg, s, e, res = 1, alpha = 8, beta = 1, thresh = 3, time = TRUE) {
  start_time <- Sys.time()
  eeg <- eeg %>%
    partition.eeg(s, e) %>%
    lower.res(res)

  analysis <- capa.mv(eeg@data[-1], type = "mean")

  canoms <- collective_anomalies(analysis) %>%
    dplyr::filter(mean.change >= alpha) %>%
    set.timevars(data = eeg@data) %>%
    format.collectives(thresh) %>%
    as_tibble()

  panoms <- point_anomalies(analysis) %>%
    dplyr::filter(strength >= beta) %>%
    set.timevars(data = eeg@data) %>%
    as_tibble()


  results <- new("analysis",
    canoms = canoms,
    panoms = panoms,
    eeg = eeg
  )
  end_time <- Sys.time()
  if (time == TRUE) {
    print(paste("Analysis completed in ", seconds_to_period(end_time - start_time)))
  }
  return(results)
}



#' Performs stepwise CAPA analysis on EEG data, saving
#' result plots for each step and writing a .csv data file
#' with every epoch/subepoch pair containing an anomaly.
#'
#' @param eeg An eeg object.
#' @param step_size Size in seconds of the sequence of data analyized
#' on each step.
#' @param res Resolution at which to perform the analyses
#' @param alpha Threshold of strength significance for collective anomalies
#' @param beta Threshold of strength significance for point anomalies
#' @param thresh How many seconds anomaly n must be from anomaly (n - 1) to
#'  consider them part of a same cluster?


analyize.stepwise <- function(eeg, step_size, res, alpha = 8, beta = 1, thresh = 3) {
  start_time <- Sys.time()
  epoch_data <- create.epoch.data()
  s <- eeg@data$Time[1]
  e <- s + step_size
  eeg_duration <- tail(eeg@data$Time, n = 1) - s
  steps <- eeg_duration %/% step_size
  if (eeg_duration %% step_size != 0) {
    steps <- steps + 1
    r <- eeg_duration %% step_size
  }

  pb <- txtProgressBar(
    min = s, max = tail(eeg@data$Time, n = 1),
    style = 3
  ) # Progress bar

  for (x in 1:(steps)) { # One step already done defining base
    if (e > tail(eeg@data, n = 1)) {
      e <- s + r
    }

    analysis <- analyze(eeg, s, e, res, alpha, beta = beta, thresh = thresh, time = FALSE)
    if (has.anomalies(analysis)) {
      plot <- plot(analysis, save = TRUE)
    }
    epoch_data <- update.epochs(epoch_data, analysis)
    s <- e
    e <- e + step_size
    setTxtProgressBar(pb, s)
  }
  write_csv(epoch_data, "results/results.csv")
  end_time <- Sys.time()
  print(paste("Stepwise analysis completed in ", seconds_to_period(end_time - start_time)))
  return(epoch_data)
}
