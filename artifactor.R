library(ggplot2)
library(reshape2)
library(anomaly)
library(patchwork)
library(dplyr)
library(cowplot)
library(gridExtra)
library(grid)
library(ggplotify)
source('dman.R')



#' Given an anomaly database, an anomaly index and a threshold in seconds,
#' returns the time at which the anomaly occurred plus the threshold in seconds.
#' The function is designed to be used in the context of removing anomalies
#' that are statistically different but clustered in a short time-span
#' (see function join_anomalous_clusters).
#' 
#' @param anom_dfr Anomaly database
#' @param index Anomaly index
#' @param threshold Threshold in seconds
#'
#' @return
thresh <- function(anom_df, index, threshold)
{
  time = anom_df$Time[index]
  time_in_secs = period_to_seconds(hms(time))
  return (time_in_secs + threshold)
}

#' Returns the subset of an anomaly database that excludes
#' all but one of the anomalies that are separated from its predecessor
#' by a time distance of x. The function aims to solve the problem of
#' massive databases resulting from CAPA analysis where a great number
#' of cases belong to a unique cluster or anomalous segment.
#' 
#' @param anom_df Anomaly database
#' @param x How many seconds from its preceding anomaly to be excluded.
#'
#' @return
join_anomalous_clusters <- function(anom_df, x, progress_bar=TRUE)
{
  df_length <- nrow(anom_df) 

    if (df_length <= 2)
  {
    return (anom_df)
    }
  
  # Progress bar
  pb <- txtProgressBar(min = 2, max = df_length, style = 3)
  
  # Vars
  anom_df <- anom_df[order(anom_df$variate),]
  keep <- c(1)
  variate <- anom_df$variate[1]
  cluster_lengths <- c()
  cluster_length <- 0
  
  for (n in 2:df_length)
  {
    if (anom_df$variate[n] != variate)
    {
      keep <- append(keep, n)
      variate <- variate + 1
      cluster_lengths <- append(cluster_lengths, cluster_length)
      cluster_length <- 0
      next
    }
    if (progress_bar == TRUE)
    {
      setTxtProgressBar(pb, n)
    }

    time <- anom_df$Time[n]
    time_in_secs <- period_to_seconds(hms(time))
    # Threshold of previous anomaly:
    threshold <- thresh(anom_df, n - 1, x) 
    if (threshold - time_in_secs >= 0 | time == anom_df$Time[n-1])
    {
      cluster_length <- cluster_length + 1
      next
    }
    cluster_lengths <- append(cluster_lengths, cluster_length)
    cluster_length <- 0
    keep <- append(keep, n)
  }
  
  if (df_length %in% keep) # If last anomaly was not in cluster
  {
    cluster_lengths <- append(cluster_lengths, 0)
  }else
  {
    cluster_lengths <- append(cluster_lengths, cluster_length)
  }
  removed <- anom_df[keep, ]
  removed$end <- anom_df$end[keep + cluster_lengths]
  
  mean_significances <- c()
  for (n in 1:nrow(removed))
  {
    mean_s <- max(anom_df$mean.change[keep[n]:(keep[n] + cluster_lengths[n])])
    mean_significances <- append(mean_significances, mean_s)
  }
  
  removed$mean.change <- mean_significances
  return (removed)
}

add_anomaly_time <- function(col_anom_df, data)
{
  
  s <- c()
  for (x in 1:nrow(col_anom_df))
  {
    start <- col_anom_df$start[x]
    time <- data$Time[start]
    s <- append(s, time)
  }
  
  s <- seconds_to_period(s)
  formated_time <- paste(round(s@hour, 4), round(minute(s), 4), round(second(s)), sep = ':')
  
  return (formated_time)
}

run_local_analysis <- function(df, start, end, resolution=1, alpha=1, thresh=3, plot=TRUE, save_image=FALSE)
{
  df <- partition_eeg_data(df, start, end)
  df <- lower_resolution(df, resolution)
  
  analysis <- capa.mv(df[-1], type='mean')
  anoms <- collective_anomalies(analysis)
  sig_anoms <- subset(anoms, anoms$mean.change >= alpha)
  sig_anoms$Time <- add_anomaly_time(sig_anoms, df)
  anoms$Time <- add_anomaly_time(anoms, df)
  sig_anoms_reduced <- join_anomalous_clusters(sig_anoms, thresh, progress_bar = FALSE)
  
  
  View(anoms)
  View(sig_anoms)
  View(sig_anoms_reduced)
  if (plot == TRUE)
  {
    point <- point_anomalies(analysis)
    anom_plot <- plot_all_channels(df, sig_anoms_reduced, point, alpha)
    print(anom_plot)
  }
  if (save_image == TRUE)
  {
    point <- point_anomalies(analysis)
    anom_plot <- plot_all_channels(df, sig_anoms_reduced, point, alpha)
    print('Saving')
    time <- paste(start, end, sep = ' to ')
    ggsave(paste('images/', time, '.png', sep = ''), plot=anom_plot)
  }
}

get_col_anoms <- function(df, resolution=1, alpha=1)
{
  if (resolution != 1)
  {
    start_time <- Sys.time() 
    df <- lower_resolution(df, resolution)
    end_time <- Sys.time()
  }

  analysis <- capa.mv(na.omit(df)[-1], type='mean') #[-1] excludes Time column.
  an_obj <- collective_anomalies(analysis)
  if (nrow(an_obj) == 0) # What if the col obj is empty (no anomalies)?
  {
    return (anl_obj)
  }
  sig_anoms <- subset(an_obj, an_obj$mean.change > alpha)
    
  if (nrow(sig_anoms) == 0) # What if the col obj is empty (no anomalies)?
  {
    return (sig_anoms)
  }
    
  sig_anoms$Time <- add_anomaly_time(sig_anoms, df)
  return (join_anomalous_clusters(sig_anoms))
}


detect_artifacts <- function(df, step_size, resolution, alpha=1) # Step size in seconds
{
  base <- get_col_anoms(partition_eeg_data(df, 1, step_size))
  View(base)
  eeg_duration <- df$Time[nrow(df)]
  s <- 1 + step_size
  e <- s + step_size - 1
  steps <- nrow(df) %/% step_size
  
  pb <- txtProgressBar(min = s, max = eeg_duration, style = 3) # Progress bar
  
  for (x in 1:steps) # One step already done defining base
  {
    #print(paste("Consider ", s, e, sep = ' '))
    if (e > df$Time[nrow(df)]) # If this is the last step...
    {
      break # For now only
      e <- nrow(df)
    }
    partitioned <- partition_eeg_data(df, s, e)
    #print(paste('Beginning analysis from ', s, ' to ', e))
    col_anoms <- get_col_anoms(partitioned, resolution, alpha)
    s <- e + 1
    e <- e + step_size
    base <- rbind(base, anoms)
    setTxtProgressBar(pb, s)
    #print(paste('Progress: ', round((s/eeg_duration) * 100, 3), '%', sep = ''))
  }
  
  return(base)
}

save_epoch_plots <- function(df, step_size, resolution, alpha=1) # Step size in seconds
{

  eeg_duration <- df$Time[nrow(df)]
  s <- 1
  e <- s + step_size - 1
  steps <- nrow(df) %/% step_size
  
  pb <- txtProgressBar(min = s, max = eeg_duration, style = 3) # Progress bar
  
  for (x in 1:steps) # One step already done defining base
  {
    #print(paste("Consider ", s, e, sep = ' '))
    if (e > df$Time[nrow(df)]) # If this is the last step...
    {
      break # For now only
      e <- nrow(df)
    }
    run_local_analysis(data, s, e, resolution, alpha, thresh = 3, plot = FALSE, save_image = TRUE)
    s <- e + 1
    e <- e + step_size
    setTxtProgressBar(pb, s)
  }
  
  return(base)
}

save_epoch_plots(data, 60, 10, 8)



plot_channel <- function(o_df, col_anoms_df, point_anoms_df, channel, alpha)
{
  
  # Sigmoid function f(alpha) to determine transparency
  # according to significance.
  thickness = round(0.5/(1+exp(-(alpha - alpha*0.25))), 3)
  
  o_df <- o_df[-1] # Exclude time variable
  col_anoms_df <- subset(col_anoms_df, col_anoms_df$variate == channel)
  point_locs <- subset(point_anoms_df, point_anoms_df$variate == channel)$location
  
  if (nrow(col_anoms_df) == 0 & length(point_locs) == 0)
  {
    plot <- ggplot(o_df, aes(1:nrow(o_df), unlist(o_df[channel]))) +
      geom_line() +
      xlab('') + ylab(channel)
  }else if (nrow(col_anoms_df) == 0)
  {
    points <- data.frame(A=point_locs, B=o_df[point_locs, channel])
    
    plot <- ggplot(o_df, aes(1:nrow(o_df), unlist(o_df[channel]))) +
      geom_line() +
      geom_point(data=points, aes(A, B), inherit.aes = FALSE, color='red') +
      xlab('') + ylab(channel)
  }else
  {
    areas <- data.frame(xmin=col_anoms_df$start, xmax=col_anoms_df$end,
                        ymin=-300, ymax=300, alpha=thickness)
    points <- data.frame(A=point_locs, B=o_df[point_locs, channel])
    
    plot <- ggplot(o_df, aes(1:nrow(o_df), unlist(o_df[channel]))) +
      geom_line() +
      geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
                alpha=thickness, fill='red',
                data=transform(areas, as.character(1:nrow(areas))),
                inherit.aes = FALSE) +
      geom_point(data=points, aes(A, B), inherit.aes = FALSE, color='red') +
      xlab('') + ylab(channel)
  }
  
  return (plot)
}

plot_all_channels <- function(o_df, col_anoms_df, point_anoms_df, alpha)
{
  plots <- list()
  for (channel in union(unique(col_anoms_df$variate), unique(point_anoms_df$variate)))
  {
    p <- plot_channel(o_df, col_anoms_df, point_anoms_df, channel, alpha)
    plots[[channel]] <- p
  }
  return(plot_grid(plotlist = plots))
}

#################

data <- na.omit(test_data)
data <- lower_resolution(data, 1)
rownames(data) <- 1:nrow(data)
View(data)


a <- detect_artifacts(data, 60, 1)
clean_a <- join_anomalous_clusters(a, 3)
View(a)
View(clean_a)



########################

run_local_analysis(data, 40, 100, resolution = 10, alpha = 8, thresh = 2, plot = FALSE, save_image = TRUE)
run_local_analysis(data, 0, 3600, resolution = 10, alpha=8, thresh=3, plot=TRUE)

run_local_analysis(data, 660, 720, resolution = 10, alpha = 8, thresh = 2, plot = TRUE, save_image = FALSE)

eleventh_minute <- partition_eeg_data(data, 660, 720)
eleventh_minute <- lower_resolution(eleventh_minute, 10)
test <- capa.mv(eleventh_minute[2:9], type='mean')
cols <- subset(collective_anomalies(test), collective_anomalies(test)$mean.change > 1)
cols$Time <- add_anomaly_time(cols, eleventh_minute)
cols_r_2 <- join_anomalous_clusters(cols, 2)
point_anoms <- point_anomalies(test)

plot_channel(eleventh_minute, cols_r, point_anoms, 1, 1) /
  plot_channel(eleventh_minute, cols_r, point_anoms, 2, 1)

fourth_minute <- partition_eeg_data(data, 240, 300)
fourth_minute <- lower_resolution(fourth_minute, 10)
test <- capa.mv(fourth_minute[2:9], type='mean')
cols <- subset(collective_anomalies(test), collective_anomalies(test)$mean.change > 1)
cols$Time <- add_anomaly_time(cols, fourth_minute)
cols_r_2 <- join_anomalous_clusters(cols, 2)
point_anoms <- point_anomalies(test)

plot_channel(eleventh_minute, cols_r_2, point_anoms, 1, 1) /
  plot_channel(eleventh_minute, cols_r_2, point_anoms, 2, 1)

plot_all_channels(eleventh_minute, cols_r_2, point_anoms, 8)
