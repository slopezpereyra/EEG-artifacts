library(ggplot2)
library(cowplot)
library(patchwork)

draw_eeg <- function(o_df, channel)
{
  o_df <- o_df[-1]
  plot <- ggplot(o_df, aes(1:nrow(o_df), unlist(o_df[channel]))) +
    geom_line()
  
  return (plot)
}

draw_col_anoms <- function(o_df, col_anoms_df, channel, alpha)
{
  
  col_anoms_df <- subset(col_anoms_df, col_anoms_df$variate == channel)
  
  # Sigmoid function f(alpha) to determine transparency
  # according to significance.
  thickness = round(0.5/(1+exp(-(alpha - alpha*0.25))), 3)
  
  areas <- data.frame(xmin=col_anoms_df$start, xmax=col_anoms_df$end,
                      ymin=-300, ymax=300, alpha=thickness)
  eeg <- draw_eeg(o_df, channel)
  plot <- eeg +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
                  alpha=thickness, fill='red',
                  data=transform(areas, as.character(1:nrow(areas))),
                  inherit.aes = FALSE) +
    xlab('') + ylab(channel)
  
  return (plot)
  
}

draw_point_anoms <- function(o_df, point_anoms_df, channel)
{ # Exclude time variable
  point_locs <- subset(point_anoms_df, point_anoms_df$variate == channel)$location
  point_values <- unlist(o_df[point_locs, channel])
  points <- tibble(A=point_locs,B=point_values)
  
  eeg <- draw_eeg(o_df, channel)
  plot <- eeg +
    geom_point(data=points, aes(A, B), inherit.aes = FALSE, color='red') +
    xlab('') + ylab(channel)
  
  return (plot)
}

draw_anomalies <- function(o_df, col_anoms_df, point_anoms_df, channel, alpha)
{
  point_locs <- subset(point_anoms_df, point_anoms_df$variate == channel)$location
  point_values <- unlist(o_df[-1][point_locs, channel])
  points <- tibble(A=point_locs,B=point_values)
  
  col_anoms_df <- subset(col_anoms_df, col_anoms_df$variate == channel)
  # Sigmoid function f(alpha) to determine transparency
  # according to significance.
  thickness = round(0.5/(1+exp(-(alpha - alpha*0.25))), 3)
  areas <- data.frame(xmin=col_anoms_df$start, xmax=col_anoms_df$end,
                      ymin=-300, ymax=300, alpha=thickness)
  
  eeg <- draw_eeg(o_df, channel)
  plot <- eeg +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
              alpha=thickness, fill='red',
              data=transform(areas, as.character(1:nrow(areas))),
              inherit.aes = FALSE) +
    geom_point(data=points, aes(A, B), inherit.aes = FALSE, color='red') +
    xlab('') + ylab(channel)
  
  return(plot)
  
}

plot_channel <- function(o_df, col_anoms_df, point_anoms_df, channel, alpha)
{
  col_anoms_df <- subset(col_anoms_df, col_anoms_df$variate == channel)
  point_locs <- subset(point_anoms_df, point_anoms_df$variate == channel)$location
  
  if (nrow(col_anoms_df) == 0 & length(point_locs) == 0)
  {
    return(draw_eeg(o_df, channel))
  }else if (nrow(col_anoms_df) == 0)
  {
    return (draw_point_anoms(o_df, point_anoms_df, channel))
  }else
  {
    return(draw_anomalies(o_df, col_anoms_df, point_anoms_df, channel, alpha))
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

