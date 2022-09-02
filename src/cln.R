# Scripts on this file are specially dedicated to cleaning
# collective outlier databases resulting from CAPA analysis.

library(lubridate)


# " Given an anomaly database, an anomaly index and a threshold in seconds,
# " returns the time at which the anomaly occurred plus the threshold in seconds.
# " The function is designed to be used in the context of removing anomalies
# " that are statistically different but clustered in a short time-span
# " (see function join_anomalous_clusters).
# "
# " @param anom_dfr Anomaly database
# " @param index Anomaly index
# " @param threshold Threshold in seconds
# "
# " @return
thresh <- function(anom_df, index, threshold) {
  time <- anom_df$Time[index]
  time_in_secs <- period_to_seconds(hms(time))
  return(time_in_secs + threshold)
}

# " Given an anomaly database, an anomaly index and a threshold in seconds,
# " returns the time at which the anomaly occurred plus the threshold in seconds.
# " The function is designed to be used in the context of removing anomalies
# " that are statistically different but clustered in a short time-span
# " (see function join_anomalous_clusters).
# "
# " @param anom_dfr Anomaly database
# " @param index Anomaly index
# " @param threshold Threshold in seconds
# "
# " @return
format.time <- function(s) {
  formated_time <- paste(round(s@hour, 4), round(minute(s), 4),
    round(second(s)),
    sep = ":"
  )
  return(formated_time)
}


# " Given a collective anomalies database, adds the formated
# " time at which each anomaly occurred based on its position
# " in the original EEG data.
# "
# " @param Collective anomalies data frame
# " @param index Data where anomalies were detected
# "
# " @return
get.time <- function(anom_df, data) {
  if (nrow(anom_df) == 0) {
    return()
  }
  s <- c()
  for (x in seq_len(nrow(anom_df))) {
    start <- anom_df[1][x, ]
    time <- data$Time[start]
    s <- append(s, time)
  }

  formated_time <- format.time(seconds_to_period(s))
  return(formated_time)
}

set_anomaly_epoch <- function(anom_df) {
  if (nrow(anom_df) == 0) {
    return(anom_df)
  }

  epoch <- anom_df$Time %>%
    hms() %>%
    period_to_seconds() %/% 30

  second <- anom_df$Time %>%
    hms() %>%
    second()
  subepoch <- ((second - 30 * (second %/% 30)) %/% 5) + 1

  anom_df$Epoch <- epoch
  anom_df$Subepoch <- subepoch
  return(anom_df)
}

set.timevars <- function(anom_df, data) {
  if (nrow(anom_df) == 0) {
    return(anom_df)
  }
  anom_df$Time <- get.time(anom_df, data)
  anom_df <- set_anomaly_epoch(anom_df)
  return(anom_df)
}

# " Returns the maximum mean change from the range of mean
# " changes that occurred in an anomalous cluster.
# "
# " @param anom_df Collective anomalies dataframe
# " @param start_points Start indexes of all clusters in anom_df
# " @param threshold End indexes of all clusters in anom_df
# "
# " @return Vector containing maximum mean change of each anomaly in anom_df
get.maxchange <- function(anom_df, start_points, end_points) {
  maxs <- list()
  for (n in seq_along(start_points)) {
    start <- start_points[n]
    end <- end_points[n]
    mean_changes <- anom_df$mean.change[start:end]
    maxs[[length(maxs) + 1]] <- round(max(mean_changes), 4)
  }
  # The values of get.maxchange are passed to cleaned
  # mean.change columns in anomaly data frames. Unlisting
  # is therefore necessary to avoid the column type is the
  # same in data frames that underwent cleaning and data frames
  # that did not (for example, those with only one anomaly).
  return(unlist(maxs))
}


# " Given the ocurrence time of an anomaly, the ocurrence time
# " of the previous anomaly and a certain threshold, determine
# " if both anomalies can be considered to be in the same cluster.
# "
# " @param previous_anom_time Formated time of (n-1)th anomaly.
# " @param time Formated time of nth anomaly.
# " @param threshold Integer determining how many seconds away two
# " anomalies must be to not be considered part of the same cluster.
# "
# " @return bool
is.clustered <- function(previous_anom_time, time, threshold) {
  time_in_secs <- period_to_seconds(hms(time))
  return(threshold - time_in_secs >= 0 | time == previous_anom_time)
}

last.cluster <- function(in_cluster, clusters, start, end) {
  if (in_cluster == TRUE) {
    clusters[[length(clusters) + 1]] <- c(start, end)
  } else {
    clusters[[length(clusters) + 1]] <- c(start, end - 1)
    clusters[[length(clusters) + 1]] <- c(end, end)
  }
  return(clusters)
}


detect.clusters <- function(anom_df, cluster_threshold) {

  # clusters[[n]][m] accesses mth element in nth cluster
  clusters <- list()
  start <- 1

  for (n in 2:nrow(anom_df)) {
    time <- anom_df$Time[n]
    predecessor_time <- anom_df$Time[n - 1]
    threshold <- thresh(anom_df, n - 1, cluster_threshold)
    in_cluster <- is.clustered(predecessor_time, time, threshold)

    if (n == nrow(anom_df)) { # Special clause needed for last element.

      clusters <- last.cluster(in_cluster, clusters, start, n)
      break
    }
    if (in_cluster == FALSE) { # If anomaly in cluster
      clusters[[length(clusters) + 1]] <- c(start, n - 1)
      start <- n
    }
  }
  return(clusters)
}

join.clusters <- function(anom_df, clusters) {
  cluster_starts <- unlist(lapply(clusters, `[[`, 1))
  cluster_ends <- unlist(lapply(clusters, `[[`, 2))

  joined <- anom_df[cluster_starts, ]
  joined$end <- anom_df$end[cluster_ends]
  joined$mean.change <- get.maxchange(anom_df, cluster_starts, cluster_ends)

  return(joined)
}

format.collectives <- function(anom_df, cluster_thresh) {
  if (nrow(anom_df) == 0) {
    return(anom_df)
  }
  anom_df <- anom_df[order(anom_df$variate), ]
  subsets <- list()
  for (channel in unique(anom_df$variate)) {
    channel_anoms <- filter(anom_df, variate == channel)
    if (nrow(channel_anoms) <= 2) {
      subsets[[length(subsets) + 1]] <- channel_anoms
      next
    }
    clusters_pos <- detect.clusters(channel_anoms, cluster_thresh)
    joined_channel_anoms <- join.clusters(channel_anoms, clusters_pos)
    subsets[[length(subsets) + 1]] <- joined_channel_anoms
  }
  return(subsets %>% reduce(full_join))
}
