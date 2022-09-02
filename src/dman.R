library(tidyverse)


# Scripts on this file are dedicated to partitioning
# and handling EEG data at the most basic level
# (sub-setting by record time, lowering data res,
# etc.).

#
partition.eeg <- function(df, start, end) {
  s_ind <- which(df$Time == start)
  e_ind <- which(df$Time == end)

  return(df[s_ind:e_ind, ])
}

lower.res <- function(df, steps) {
  i <- c(1) # Include one since following computation skips it.
  if (steps == 1) {
    return(df)
  }
  for (n in seq_len(nrow(df) / steps)) {
    i <- append(i, n * steps + 1)
  }

  return(df[i, ])
}

create.epoch.data <- function() {
  results <- tibble(
    epoch = numeric(),
    subepoch = numeric(),
    channels = list(),
    segment_strength = numeric(),
    point_strength = numeric()
  )
  return(results)
}

update.epochs <- function(epoch_data, analysis) {
  colls <- analysis[[1]]
  points <- analysis[[2]]

  collective_epochs <- unique(paste(colls$Epoch, colls$Subepoch))
  point_epochs <- unique(paste(points$Epoch, points$Subepoch))

  epochs <- union(collective_epochs, point_epochs)
  epoch_list <- lapply(strsplit(epochs, " "), as.numeric)
  chans <- union(unique(colls$variate), unique(points$variate))

  epoch_data <- add_row(epoch_data,
    epoch = unlist(lapply(epoch_list, `[[`, 1)),
    subepoch = unlist(lapply(epoch_list, `[[`, 2))
  )

  return(epoch_data)
}
