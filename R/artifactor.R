# Scripts on this file are dedicated to performing M-CAPA analysis
# on EEG data.

library(anomaly)
library(methods)
source("R/cln.R")
source("R/analysis.r")
source("R/eeg.R")


# " Perform M-CAPA analysis on EEG data in a given range of seconds
# " and return results dplyr::filtered by anomaly strength.
# "
# " @param df EEG data
# " @param s First second of timespan to analyze
# " @param e Last second of timespan to analyze
# " @param res Resolution at which to perform analysis
# "            (see lower.resolution function)
# " @param alpha Threshold of strength significance for collective anomalies
# " @param beta Threshold of strength significance for point anomalies
# " @param thresh How many seconds anomaly n must be from anomaly (n - 1) to
# "               consider them part of a same cluster?
# " @param save_origin Return df along with results?

# " @return List containing collective anomaly data, point anomaly data
# " @       and df if save_origin is true.
analyze <- function(eeg, s, e, res = 1, alpha = 1, beta = 1, thresh = 3) {
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
    origin = eeg@data
  )
  return(results)
}



# " Performs stepwise M-CAPA analysis on EEG data, saving
# " result plots for each step and writes a .csv data file
# " with every epoch/subepoch pair containing an anomaly.
# "
# " @param df EEG data
# " @param s First second of timespan to analyze
# " @param e Last second of timespan to analyze
# " @param res Resolution at which to perform analysis
# "            (see lower.resolution function)
# " @param alpha Threshold of strength significance for collective anomalies
# " @param beta Threshold of strength significance for point anomalies
# " @param thresh How many seconds anomaly n must be from anomaly (n - 1) to
# "               consider them part of a same cluster?
# " @param save_origin Return df along with results?

# " @return List containing collective anomaly data, point anomaly data
# " @       and df if save_origin is true.

analyize.stepwise <- function(eeg, step_size, res, alpha = 1, beta = 1) {
  eeg_duration <- tail(eeg@data$Time, n = 1)
  s <- 1
  e <- s + step_size - 1
  steps <- nrow(eeg@data) %/% step_size
  epoch_data <- create.epoch.data()

  pb <- txtProgressBar(min = 1, max = eeg_duration, style = 3) # Progress bar

  for (x in 1:steps) { # One step already done defining base
    if (e > tail(eeg@data$Time, n = 1)) { # If this is the last step...
      break # For now only
      e <- nrow(df)
    }
    analysis <- analyze(eeg, s, e, res, alpha, beta = beta, thresh = 3)
    if (has.anomalies(analysis)) {
      plot <- plot(analysis, save = TRUE)
    }
    epoch_data <- update.epochs(epoch_data, analysis)
    s <- e
    e <- e + step_size
    setTxtProgressBar(pb, s)
  }
  write_csv(epoch_data, "results.csv")
}
