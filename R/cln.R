# Scripts on this file are specially dedicated to cleaning
# collective outlier data frames resulting from CAPA analysis.
# (Functions are not dedicated to point anomalies data frames).

library(lubridate)

# " Given an collective anomalies (or canoms) data frame,
# " an index  i and a threshold in seconds, returns the time at which
# ' anomaly canoms[i] occurred plus the threshold in seconds.
# "
# " @param canoms Collective anomalies data frame, such as the @canoms attribute in an analysis object.
# " @param i Anomaly index
# " @param threshold Threshold in seconds
# "
# " @return An integer. The time of occurence of anomaly canoms[i] plus threshold.
thresh <- function(canoms, i, threshold) {
  time <- canoms$Time[i]
  time_in_secs <- period_to_seconds(time)
  return(time_in_secs + threshold)
}

# " Given a lubridate time object s, returns a character formatted in the
# " hours:minutes:seconds.
# "
# " @param s A lubridate time object.
# "
# " @return A character.
format.time <- function(s) {
  formated_time <- paste(round(s@hour, 4), round(minute(s), 4),
    round(second(s)),
    sep = ":"
  )
  return(formated_time)
}


# " Given a canoms or panoms data frame, returns a vector containing
# " the lubridate period objects corresponding to the time of occurrence,
# " in the EEG record, of every anomaly of the data frame.
# "
# " @param anoms A canoms or panoms data frame.
# "
# " @return A vector of lubridate periods
get.time <- function(anoms, data) {
  return(seconds_to_period(data$Time[unlist(anoms[1])]))
}

# " Given a canoms or panoms data frame, returns an identical
# " canoms or panoms with the inclusion of epoch and subepoch columns.
# "
# " @param anoms A canoms or panoms data frame.
# "
# " @return A data frame
set_anomaly_epoch <- function(anoms) {
  if (nrow(anoms) == 0) {
    return(anoms)
  }

  epoch <- anoms$Time %>%
    period_to_seconds() %/% 30

  second <- anoms$Time %>%
    second()
  subepoch <- ((second - 30 * (second %/% 30)) %/% 5) + 1

  anoms$Epoch <- epoch
  anoms$Subepoch <- subepoch
  return(anoms)
}

# " Given a canoms or panoms data frame and its origin, returns an identical
# " canoms or panoms data frame with the inclusion of all time variables:
# " Time, Epoch and Subepoch.
# "
# " @param anoms A canoms or panoms data frame.
# " @data The data where anomalies were detected.
# " @return A data frame
set.timevars <- function(anoms, data) {
  if (nrow(anoms) == 0) {
    return(anoms)
  }
  anoms$Time <- get.time(anoms, data)
  anoms <- set_anomaly_epoch(anoms)
  return(anoms)
}

# " Returns the maximum of all mean changes that occurred in
# " an anomalous cluster.
# "
# " @param anom_df A canoms dataframe.
# " @param start_points Start indexes of all clusters in canoms
# " @param threshold End indexes of all clusters in anoms
# "
# " @return Vector containing maximum mean change of each anomaly in canoms
get.maxchange <- function(canoms, start_points, end_points) {
  maxs <- list()
  for (n in seq_along(start_points)) {
    start <- start_points[n]
    end <- end_points[n]
    mean_changes <- canoms$mean.change[start:end]
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
  time_in_secs <- period_to_seconds(time)
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

# " Given a single-channel canoms data frame, returns a list of all
# "  clusters of anomalies. Each element in the list is a pair of
# " the form (s, e) containing the start and the end indexes of each cluster.
# "
# " @param canoms A canoms data frame with a single channel.
# "
# " @return A list of (s, e) pairs.
detect.clusters <- function(canoms, cluster_threshold) {
  clusters <- list()
  start <- 1

  for (n in 2:nrow(canoms)) {
    time <- canoms$Time[n]
    predecessor_time <- canoms$Time[n - 1]
    threshold <- thresh(canoms, n - 1, cluster_threshold)
    in_cluster <- is.clustered(predecessor_time, time, threshold)

    if (n == nrow(canoms)) { # Special clause needed for last element.

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

# " Given a single-variate canoms data frame, returns an equal size
# " or smaller data frame where all anomalies belonging to the same
# " cluster were reduced to a single anomaly.
# "
# " @param canoms A canoms data frame with a single channel.
# " @param clusters A list of (s, e) pairs containing start and end indexes of all clusters in canoms.
# "
# " @return A canoms data frame.
join.clusters <- function(canoms, clusters) {
  cluster_starts <- unlist(lapply(clusters, `[[`, 1))
  cluster_ends <- unlist(lapply(clusters, `[[`, 2))

  joined <- canoms[cluster_starts, ]
  joined$end <- canoms$end[cluster_ends]
  joined$mean.change <- get.maxchange(canoms, cluster_starts, cluster_ends)

  return(joined)
}

# " Given a multi-channel canoms data frame, returns a formatted version of the data frame
# " where all anomalies belonging to the same clusterwere reduced to a single
# " anomaly in every record channel.
# "
# " @param canoms A canoms data frame.
# " @param clusters A list of (s, e) pairs containing start and end indexes of all clusters in canoms.
# "
# " @return A canoms data frame.
format.collectives <- function(anom_df, cluster_thresh) {
  if (nrow(anom_df) == 0) {
    return(anom_df)
  }
  anom_df <- anom_df[order(anom_df$variate), ]
  subsets <- list()
  for (channel in unique(anom_df$variate)) {
    channel_anoms <- dplyr::filter(anom_df, variate == channel)
    if (nrow(channel_anoms) <= 2) {
      subsets[[length(subsets) + 1]] <- channel_anoms
      next
    }
    clusters_pos <- detect.clusters(channel_anoms, cluster_thresh)
    joined_channel_anoms <- join.clusters(channel_anoms, clusters_pos)
    subsets[[length(subsets) + 1]] <- joined_channel_anoms
  }
  suppressMessages()
  return(subsets %>% reduce(full_join))
}
