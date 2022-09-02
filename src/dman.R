# Scripts on this file are dedicated to partitioning
# and handling EEG data at the most basic level
# (sub-setting by record time, lowering data res,
# etc.).


library(tidyverse)

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
  canoms <- analysis@canoms
  panoms <- analysis@panoms

  collective_epochs <- unique(paste(canoms$Epoch, canoms$Subepoch))
  point_epochs <- unique(paste(panoms$Epoch, panoms$Subepoch))

  epochs <- union(collective_epochs, point_epochs)
  epoch_list <- lapply(strsplit(epochs, " "), as.numeric)
  chans <- union(unique(canoms$variate), unique(panoms$variate))

  epoch_data <- add_row(epoch_data,
    epoch = unlist(lapply(epoch_list, `[[`, 1)),
    subepoch = unlist(lapply(epoch_list, `[[`, 2))
  )

  return(epoch_data)
}
