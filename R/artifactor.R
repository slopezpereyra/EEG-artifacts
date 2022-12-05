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
#' @param s int First second of timespan to analyze
#' @param e int Last second of timespan to analyze
#' @param res int Resolution at which to perform analysis
#' @param alpha float Threshold of strength significance for collective anomalies
#' @param beta float Threshold of strength significance for point anomalies
#' @param time bool Show process time?
#' @return An analysis object
#' @export
#'
analyze <- function(eeg, s = -1, e = -1, res = 1, alpha = 8, beta = 1, time = TRUE) {
    s <- ifelse(s != -1, s, eeg@data$Time[1])
    e <- ifelse(e != -1, e, tail(eeg@data$Time, n = 1))
    start_time <- Sys.time()
    eeg <- eeg %>%
        subset_eeg(s, e) %>%
        resample_eeg(res)

    analysis <- capa.mv(eeg@data[-1], type = "mean")

    canoms <- collective_anomalies(analysis) %>%
        dplyr::filter(mean.change >= alpha) %>%
        set_timevars(data = eeg@data) %>%
        as_tibble()

    panoms <- point_anomalies(analysis) %>%
        dplyr::filter(strength >= beta) %>%
        set_timevars(data = eeg@data) %>%
        as_tibble()


    results <- new("analysis",
        canoms = canoms,
        panoms = panoms,
        eeg = eeg
    )
    end_time <- Sys.time()
    if (time) {
        print(paste("Analysis completed in ", (end_time - start_time), " seconds"))
    }
    return(results)
}

#' Performs stepwise (or epoch by epoch) CAPA analysis on EEG data.
#'
#' @param eeg An eeg object.
#' @param step_size Size in seconds of each step.
#' @param res Resolution at which to perform the analyses. Defaults to 1.
#' @param alpha Threshold of strength significance for collective anomalies.
#' @param beta Threshold of strength significance for point anomalies.
#' @return An analysis object.
#' @export
stepwise_analysis <- function(eeg, step_size, res = 1, alpha = 8, beta = 1, write = FALSE) {

    # Variables for stepwise Iteration.
    s <- eeg@data$Time[1]
    e <- s + step_size
    measures_per_step <- which(eeg@data$Time == s + step_size)
    eeg_duration <- tail(eeg@data$Time, n = 1) - s
    steps <- eeg_duration %/% step_size
    if (eeg_duration %% step_size != 0) {
        steps <- steps + 1
        r <- eeg_duration %% step_size
    }

    an <- analyze(eeg, s, e, alpha = alpha, beta = beta)
    s <- e
    e <- e + step_size
    for (x in 2:steps) {
        if (e > tail(eeg@data["Time"], n = 1)) {
            e <- s + r
        }
        analysis <- analyze(eeg, s, e, res, alpha, beta = beta)
        # Ensure point locations match the original EEG and no the
        # step subset.
        displacement <- (x - 1) * measures_per_step
        analysis@panoms[1] <- analysis@panoms[1] + displacement
        analysis@canoms[1] <- analysis@canoms[1] + displacement
        analysis@canoms[2] <- analysis@canoms[2] + displacement
        an <- merge(an, analysis)
        s <- e
        e <- e + step_size
    }
    return(an)
}
