library(ggplot2)
library(cowplot)
library(patchwork)

draw.egg <- function(o_df, channel) {
  o_df <- o_df[-1]
  plot <- ggplot(o_df, aes(1:nrow(o_df), unlist(o_df[channel]))) +
    geom_line()

  return(plot)
}

draw.clusters <- function(o_df, col_anoms, channel) {
  col_anoms <- subset(col_anoms, col_anoms$variate == channel)
  areas <- data.frame(
    xmin = col_anoms$start, xmax = col_anoms$end,
    ymin = -300, ymax = 300, alpha = 0.3
  )
  eeg <- draw.egg(o_df, channel)
  plot <- eeg +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      alpha = 0.3, fill = "red",
      data = transform(areas, as.character(1:nrow(areas))),
      inherit.aes = FALSE
    ) +
    xlab("") + ylab(channel)

  return(plot)
}

draw.points <- function(o_df, point_anoms, channel) {
  point_locs <- point_anoms %>%
    filter(variate == channel) %>%
    .$location
  point_values <- unlist(o_df[-1][point_locs, channel])
  points <- tibble(A = point_locs, B = point_values)

  eeg <- draw.egg(o_df, channel)
  plot <- eeg +
    geom_point(data = points, aes(A, B), inherit.aes = FALSE, color = "red") +
    xlab("") + ylab(channel)

  return(plot)
}

draw.anomalies <- function(o_df, col_anoms, point_anoms, channel) {
  point_locs <- point_anoms %>%
    filter(variate == channel) %>%
    .$location
  point_values <- unlist(o_df[-1][point_locs, channel])
  points <- tibble(A = point_locs, B = point_values)

  col_anoms <- subset(col_anoms, col_anoms$variate == channel)
  areas <- data.frame(
    xmin = col_anoms$start, xmax = col_anoms$end,
    ymin = -300, ymax = 300, alpha = 0.3
  )

  eeg <- draw.egg(o_df, channel)
  plot <- eeg +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      alpha = 0.3, fill = "red",
      data = transform(areas, as.character(1:nrow(areas))),
      inherit.aes = FALSE
    ) +
    geom_point(data = points, aes(A, B), inherit.aes = FALSE, color = "red") +
    xlab("") + ylab(channel)

  return(plot)
}

plot.channel <- function(o_df, col_anoms, point_anoms, channel) {
  cols <- col_anoms %>% filter(variate == channel)
  point_locs <- point_anoms %>%
    filter(variate == channel) %>%
    .$location
  if (nrow(cols) == 0 & length(point_locs) == 0) {
    return(draw.egg(o_df, channel))
  } else if (nrow(cols) == 0) {
    return(draw.points(o_df, point_anoms, channel))
  } else {
    return(draw.anomalies(o_df, cols, point_anoms, channel))
  }
  return(plot)
}

plot.channels <- function(o_df, col_anoms, point_anoms) {
  plots <- list()
  channels <- union(unique(col_anoms$variate), unique(point_anoms$variate))
  for (channel in channels)
  {
    p <- plot.channel(o_df, col_anoms, point_anoms, channel)
    plots[[channel]] <- p
  }
  return(plot_grid(plotlist = plots))
}

plot.analysis <- function(analysis, save_plot = FALSE) {
  anom_plot <- plot.channels(analysis[[3]], analysis[[1]], analysis[[2]])
  start <- seconds_to_period(analysis[[3]]$Time[1]) %>% format.time()
  end <- seconds_to_period(tail(analysis[[3]]$Time, 1)) %>% format.time()
  if (save_plot == TRUE) {
    time <- paste(start, end, sep = " to ")
    ggsave(paste("/home/santi/work/EEG-artifacts/images/", time, ".png", sep = ""), plot = anom_plot)
  }
  return(anom_plot)
}
