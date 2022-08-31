library(anomaly)
source("dman.R")
source("cln.R")
source("plotter.R")

data <- read_csv("./data/test_data.txt")
data <- na.omit(data)
data[10:14] <- NULL


analyze <- function(df, start, end, res = 1, alpha = 1, thresh = 3,
                                    plot = TRUE, save_plot = FALSE) {

  df <- df %>% partition_eeg_data(start, end) %>% lower_res(res)

  analysis <- capa.mv(df[-1], type = "mean")
  sig_anoms <- collective_anomalies(analysis) %>% filter(mean.change >= alpha)
  #sig_anoms <- subset(anoms, anoms$mean.change >= alpha)

  if (nrow(sig_anoms) == 0) {
    return(sig_anoms)
  }
  sig_anoms$Time <- add_anomaly_time(sig_anoms, df)

  sig_anoms_reduced <- clean_anomaly_data(sig_anoms, thresh)

  # NEEDS REFACTORING

  if (plot == TRUE) {
    point <- point_anomalies(analysis)
    print("Plotting")
    anom_plot <- plot_all_channels(df, sig_anoms_reduced, point, alpha)
    print(anom_plot)
  }
  if (save_plot == TRUE) {
    point <- point_anomalies(analysis)
    anom_plot <- plot_all_channels(df, sig_anoms_reduced, point, alpha)
    time <- paste(start, end, sep = " to ")
    ggsave(paste("./images/", time, ".png", sep = ""), plot = anom_plot)
  }
return(sig_anoms_reduced)
}


get_col_anoms <- function(df, res = 1, alpha = 1) {
  if (res != 1) {
    df <- lower_res(df, res)
  }

  analysis <- capa.mv(na.omit(df)[-1], type = "mean") #[-1] excludes Time var.
  an_obj <- collective_anomalies(analysis)
  if (nrow(an_obj) == 0) { # What if the col obj is empty (no anomalies)?
    return(an_obj)
  }
  sig_anoms <- subset(an_obj, an_obj$mean.change > alpha)

  if (nrow(sig_anoms) == 0) {# What if the col obj is empty (no anomalies)?
    return(sig_anoms)
  }

  sig_anoms$Time <- add_anomaly_time(sig_anoms, df)
  return(join_anomalous_clusters(sig_anoms))
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
