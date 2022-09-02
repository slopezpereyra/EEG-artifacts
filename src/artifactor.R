# Scripts on this file are dedicated to performing M-CAPA analysis
# on EEG data.

library(anomaly)
library(methods)
source("dman.R")
source("cln.R")
source("plotter.R")

setClass("analysis",
  slots = list(
    canoms = "data.frame",
    panoms = "data.frame",
    origin = "data.frame"
  )
)

setMethod(
  "show",
  "analysis",
  function(object) {
    View(object@canoms)
    View(object@panoms)
    View(object@origin)
  }
)


# " Check if analysis as returned by the analyze() function
# " found anomalies or not.
# " @param analysis An analysis as returned by the analyze() function.

setMethod(
  "has.anomalies",
  "analysis",
  function(object) {
    return(nrow(object@canoms) > 0 | nrow(object@panoms) > 0)
  }
)

# " Perform M-CAPA analysis on EEG data in a given range of seconds
# " and return results filtered by anomaly strength.
# "
# " @param df EEG data
# " @param start First second of timespan to analyze
# " @param end Last second of timespan to analyze
# " @param res Resolution at which to perform analysis
# "            (see lower.resolution function)
# " @param alpha Threshold of strength significance for collective anomalies
# " @param beta Threshold of strength significance for point anomalies
# " @param thresh How many seconds anomaly n must be from anomaly (n - 1) to
# "               consider them part of a same cluster?
# " @param save_origin Return df along with results?

# " @return List containing collective anomaly data, point anomaly data
# " @       and df if save_origin is true.
analyze <- function(df, start, end, res = 1, alpha = 1, beta = 1, thresh = 3) {
  df <- df %>%
    partition.eeg(start, end) %>%
    lower.res(res)
  analysis <- capa.mv(df[-1], type = "mean")

  canoms <- collective_anomalies(analysis) %>%
    filter(mean.change >= alpha) %>%
    set.timevars(data = df) %>%
    format.collectives(thresh) %>%
    as_tibble()

  panoms <- point_anomalies(analysis) %>%
    filter(strength >= beta) %>%
    set.timevars(data = df) %>%
    as_tibble()


  results <- new("analysis", canoms = canoms, panoms = panoms, origin = df)
  return(results)
}



# " Performs stepwise M-CAPA analysis on EEG data, saving
# " result plots for each step and writes a .csv data file
# " with every epoch/subepoch pair containing an anomaly.
# "
# " @param df EEG data
# " @param start First second of timespan to analyze
# " @param end Last second of timespan to analyze
# " @param res Resolution at which to perform analysis
# "            (see lower.resolution function)
# " @param alpha Threshold of strength significance for collective anomalies
# " @param beta Threshold of strength significance for point anomalies
# " @param thresh How many seconds anomaly n must be from anomaly (n - 1) to
# "               consider them part of a same cluster?
# " @param save_origin Return df along with results?

# " @return List containing collective anomaly data, point anomaly data
# " @       and df if save_origin is true.

analyize.stepwise <- function(df, step_size, res, alpha = 1, beta = 1) {
  eeg_duration <- df$Time[nrow(df)]
  s <- 1
  e <- s + step_size - 1
  steps <- nrow(df) %/% step_size
  epoch_data <- create.epoch.data()

  pb <- txtProgressBar(min = 1, max = eeg_duration, style = 3) # Progress bar

  for (x in 1:steps) { # One step already done defining base
    if (e > df$Time[nrow(df)]) { # If this is the last step...
      break # For now only
      e <- nrow(df)
    }

    analysis <- analyze(df, s, e, res, alpha, beta = beta, thresh = 3)
    if (has.anomalies(analysis)) {
      plot <- plot.analysis(analysis, save_plot = TRUE)
    }
    epoch_data <- update.epochs(epoch_data, analysis)
    s <- e
    e <- e + step_size
    setTxtProgressBar(pb, s)
  }
  write_csv(epoch_data, "results.csv")
}
