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
#'               collective anomaly (n - 1) to consider them part of a same cluster?

#' @return An analysis object
analyze <- function(eeg, s, e, res = 1, alpha = 8, beta = 1, thresh = 3) {
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
  eeg_duration <- tail(eeg@data$Time, n = 1)
  s <- eeg@data$Time[1]
  e <- s + step_size - 1
  steps <- nrow(eeg@data) %/% step_size
  epoch_data <- create.epoch.data()

  pb <- txtProgressBar(min = s, max = eeg_duration, style = 3) # Progress bar

  for (x in 1:steps) { # One step already done defining base
    if (e > tail(eeg@data$Time, n = 1)) { # If this is the last step...
      break # For now only
      e <- nrow(df)
    }
    analysis <- analyze(eeg, s, e, res, alpha, beta = beta, thresh = thresh)
    if (has.anomalies(analysis)) {
      plot <- plot(analysis, save = TRUE)
    }
    epoch_data <- update.epochs(epoch_data, analysis)
    s <- e
    e <- e + step_size
    setTxtProgressBar(pb, s)
  }
  write_csv(epoch_data, "results/results.csv")
}
