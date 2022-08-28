library(anomaly)
source('dman.R')
source('cln.R')
source('plotter.R')

data <- read_csv('./data/test_data.txt')
data <- na.omit(data)

run_local_analysis <- function(df, start, end, resolution=1, alpha=1, thresh=3, plot=TRUE, save_image=FALSE)
{
  df <- partition_eeg_data(df, start, end)
  df <- lower_resolution(df, resolution)
  
  analysis <- capa.mv(df[-1], type='mean')
  anoms <- collective_anomalies(analysis)
  sig_anoms <- subset(anoms, anoms$mean.change >= alpha)
  sig_anoms$Time <- add_anomaly_time(sig_anoms, df)
  anoms$Time <- add_anomaly_time(anoms, df)
  sig_anoms_reduced <- clean_anomaly_data(sig_anoms, thresh)
  
  
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

run_local_analysis(data, 660, 720, resolution = 10, alpha=8, thresh=2, plot=TRUE)

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
    return (an_obj)
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
eleventh_minute <- lower_resolution(eleventh_minute, 1)
test <- capa.mv(eleventh_minute[2:9], type='mean')
cols <- subset(collective_anomalies(test), collective_anomalies(test)$mean.change > 8)
cols$Time <- add_anomaly_time(cols, eleventh_minute)
cols <- join_anomalous_clusters(cols, 2)
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
