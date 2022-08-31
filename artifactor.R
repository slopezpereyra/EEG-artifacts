library(anomaly)
source("dman.R")
source("cln.R")
source("plotter.R")

data <- read_csv("./data/test_data.txt")
data <- na.omit(data)
data[10:14] <- NULL

#" Perform CAPA analysis on EEG data in a given range of seconds
#" and return results filtered by anomaly strength.
#"
#" @param df EEG data
#" @param start First second of timespan to analyze
#" @param end Last second of timespan to analyze
#" @param res Resolution at which to perform analysis 
#"            (see lower_resolution function)
#" @param alpha Threshold of strength significance for collective anomalies 
#" @param beta Threshold of strength significance for point anomalies
#" @param thresh How many seconds anomaly n must be from anomaly (n - 1) to
#"               consider them part of a same cluster?
#" @param end save_origin Return df along with results?

#" @return List containing collective anomaly data, point anomaly data
#" @       and df if save_origin is true.
analyze <- function(df, start, end, res = 1, alpha = 1, beta = 1, thresh = 3,
                    save_origin = TRUE) {

  df <- df %>% partition_eeg_data(start, end) %>% lower_res(res)
  analysis <- capa.mv(df[-1], type = "mean")
  col_anoms <- collective_anomalies(analysis) %>% filter(mean.change >= alpha)
  point_anoms <- point_anomalies(analysis) %>% filter(strength >= beta) 

  if (nrow(col_anoms) != 0) {
    col_anoms$Time <- add_anomaly_time(col_anoms, df)
    col_anoms <- col_anoms %>% format_collective_data(thresh)
  }

  results <- list("Collective anomalies" = col_anoms,
                  "Point anomalies" = point_anoms,
                  "Origin" = df)
  if (save_origin == FALSE) {
    results <- head(results, -1)
  }
  return(results)
}

has_anomalies <- function(analysis){
  # We can't use the any() function because analysis may contain
  # original data (non-empty always).
  return(nrow(analysis[[1]]) > 0 | nrow(analysis[[2]]) > 0)
}

view_analysis <- function(analysis){
  lapply(analysis[1:2], View)
}

analyze_by_steps <- function(df, step_size, res, alpha = 1) {
  base <- get_col_anoms(partition_eeg_data(df, 1, step_size))
  View(base)
  eeg_duration <- df$Time[nrow(df)]
  s <- 1 + step_size
  e <- s + step_size - 1
  steps <- nrow(df) %/% step_size

  pb <- txtProgressBar(min = s, max = eeg_duration, style = 3) # Progress bar

  for (x in 1:steps) {# One step already done defining base
    if (e > df$Time[nrow(df)]) { # If this is the last step...
      break # For now only
      e <- nrow(df)
    }
    partitioned <- partition_eeg_data(df, s, e)
    col_anoms <- get_col_anoms(partitioned, res, alpha)
    s <- e + 1
    e <- e + step_size
    base <- rbind(base, col_anoms)
    setTxtProgressBar(pb, s)
  }

  return(base)
}

save_epoch_plots <- function(df, step_size, res, alpha = 1) {

  eeg_duration <- df$Time[nrow(df)]
  s <- 1
  e <- s + step_size - 1
  steps <- nrow(df) %/% step_size

  pb <- txtProgressBar(min = s, max = eeg_duration, style = 3) # Progress bar

  for (x in 1:steps){ # One step already done defining base
    if (e > df$Time[nrow(df)]) { # If this is the last step...
      break # For now only
      e <- nrow(df)
    }
    analyze(data, s, e, res, alpha, thresh = 3, FALSE, TRUE)
    s <- e + 1
    e <- e + step_size
    setTxtProgressBar(pb, s)
  }
}

