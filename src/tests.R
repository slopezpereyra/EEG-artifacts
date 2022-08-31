source('artifactor.R')

data <- read_csv("./data/test_data.txt")
data <- na.omit(data)
data[10:14] <- NULL

first_half_hour <- partition_eeg_data(data, 1, 1800)
save_epoch_plots(data, 30, 1, 8)
analysis_7 <- analyze(data, 901, 960, res = 1, alpha = 8)
view_analysis(analysis_7)

half_hour_analysis <- analyze(first_half_hour, 1, 1800, res = 1, alpha = 8)
plot_analysis(half_hour_analysis)

# Step by step

save_epoch_plots(first_half_hour, 60, 1, 8)
view_analysis(inspect)
plot_analysis(inspect)
# Error in minute 4: points plotted wrongly channel 3


































#################




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
