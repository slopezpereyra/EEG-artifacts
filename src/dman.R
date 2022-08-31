# Scripts on this file are dedicated to partitioning
# and handling EEG data at the most basic level
# (sub-setting by record time, lowering data res,
# etc.).

#
partition_eeg_data <- function(df, start, end) {
  s_ind <- which(df$Time == start)
  e_ind <- which(df$Time == end)

  return(df[s_ind:e_ind, ])
}

lower_res <- function(df, steps) {
  i <- c(1) # Include one since following computation skips it.
  if (steps == 1) {
    return(df)
  }
  for (n in seq_len(nrow(df) / steps)) {
    i <- append(i, n * steps + 1)
  }

  return(df[i, ])
}