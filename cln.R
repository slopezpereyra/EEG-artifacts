# Scripts on this file are specially dedicated to cleaning
# collective outlier databases resulting from CAPA analysis.

library(tidyverse)
library(lubridate)


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


get_max_changes <- function(anom_df, start_points, end_points)
{
  maxs <- list()
  for (n in 1:length(start_points))
  {
    start <- start_points[n]
    end <- end_points[n]
    mean_changes <- anom_df$mean.change[start:end]
    maxs[[length(maxs) + 1]] <- round(max(mean_changes), 4)
  }
  # The values of get_max_changes are passed to cleaned
  # mean.change columns in anomaly data frames. Unlisting
  # is therefore necessary to avoid the column type is the
  # same in data frames that underwent cleaning and data frames
  # that did not (for example, those with only one anomaly).
  return (unlist(maxs))
}

is_in_cluster <- function(previous_anom_time, time, threshold)
{
  time_in_secs <- period_to_seconds(hms(time))
  return (threshold - time_in_secs >= 0 | time == previous_anom_time)
}

last_anomaly_cluster <- function(in_cluster, clusters, start, end)
{
  if (in_cluster == TRUE)
  {
    clusters[[length(clusters) + 1]] <- c(start, end)
  }else
  {
    clusters[[length(clusters) + 1]] <- c(start, end - 1)
    clusters[[length(clusters) + 1]] <- c(end, end)
  }
  return (clusters)
}

detect_clusters <- function(anom_df, cluster_threshold)
{
  
  # clusters[[n]][m] accesses mth element in nth cluster
  clusters <- list()
  start <- 1
  
  for (n in 2:nrow(anom_df))
  {
    time <- anom_df$Time[n]
    predecessor_time <- anom_df$Time[n-1]
    threshold <- thresh(anom_df, n - 1, cluster_threshold)
    in_cluster <- is_in_cluster(predecessor_time, time, threshold)
    
    if (n == nrow(anom_df)) # Special clause needed for last element.
    {
      clusters <- last_anomaly_cluster(in_cluster, clusters, start, n)
      break
    }
    if (in_cluster == FALSE) # If anomaly in cluster
    {
      clusters[[length(clusters) + 1]] <- c(start, n - 1)
      start <- n
    }
  }
  return (clusters)
}


join_clusters <- function(anom_df, clusters)
{
  cluster_starts <- unlist(lapply(clusters, `[[`, 1))
  cluster_ends <- unlist(lapply(clusters, `[[`, 2))
  
  joined <- anom_df[cluster_starts, ]
  joined$end <- anom_df$end[cluster_ends]
  joined$mean.change <- get_max_changes(anom_df, cluster_starts, cluster_ends)
  
  return (joined) 
}

clean_anomaly_data <- function(anom_df, cluster_thresh)
{
  anom_df <- anom_df[order(anom_df$variate), ]
  subsets <- list()
  for (channel in unique(anom_df$variate))
  {
    channel_anoms <- filter(anom_df, variate==channel)
    if (nrow(channel_anoms) <= 2)
    {
      subsets[[length(subsets)+ 1]] <- channel_anoms
      next
    }
    clusters_pos <- detect_clusters(channel_anoms, cluster_thresh)
    joined_channel_anoms <- join_clusters(channel_anoms, clusters_pos)
    subsets[[length(subsets)+ 1]] <- joined_channel_anoms
  }
  return (subsets %>% reduce(full_join))
}
