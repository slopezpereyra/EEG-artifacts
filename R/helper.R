
#' @description 
#' Reads a .edf file and returns a tibble with its EEG data.\r
#' *Important* Only columns whose name starts with EEG or EOG are kept.
#' This means the function removes snoring or EMG signals if possible.
#' @param file Path to the EDF file.
#'
#' @return tibble
#' @export
read_edf <- function(file){
    eeg <- edf::read.edf(file)
    x <- tibble::as_tibble(eeg[[3]], .name_repair = "minimal")
    t <- x[, 1][[1]][2] # Take the list of time values from first signal
    x <- x[-2, ]        # Remove the second row whose values are lists
                        # of time values (now only signals exist in x).

    # Select EEG signals and unnest the lists of the single row of x.
    # Add time colum at cool index 1.
    x  <- x %>% select(starts_with("EEG") | starts_with("EOG")) %>%
                unnest(cols = everything()) %>%
                mutate(Time = unlist(t)) %>%
                relocate(Time, .before = 1)

    return(x)
}

#' @export
znormalization <- function(x) {
    if (is.null(x) || length(x) == 0) {
        return(x)
    }
    return((x - mean(x)) / sd(x))
}

#' @export
minmax_normalization <- function(x) {
    if (is.null(x) || length(x) == 0) {
        return(x)
    }
    normalized <- (x - min(x)) / (max(x) - min(x))
    return(normalized)
}

#' @export
canoms_avg_epoch_strength <- function(df) {
    a <- df %>%
            dplyr::group_by(Epoch, Subepoch) %>%
            dplyr::summarise_at(vars(mean.change), list(cAnomStrength = mean))
    print(a)
    return(a)
}

#' Helper function
#' @export
panoms_avg_epoch_strength <- function(df) {
    a <- df %>%
            dplyr::group_by(Epoch, Subepoch) %>%
            dplyr::summarise_at(vars(strength), list(pAnomStrength = mean))
    return(a)
}
# Helper function
#' @export
set_epochs <- function(df, epoch = 30, subepochs = FALSE) {
    # Get quotient and remainder of euclidean division
    # of each time in seconds by 30.
    q <- df$Time %/% epoch
    df <- df %>% tibble::add_column(Epoch = as.factor(q + 1), .after = 1)
    if (subepochs) {
        r <- df$Time %% epoch
        qprime <- r %/% 5
        df <- df %>% tibble::add_column(Subepoch = as.factor(qprime + 1), .after = 2)
    }
    return(df)
}

#' @export
window_split <- function(x, samples_per_split){
  split(x, ceiling(seq_along(x)/samples_per_split))
}

#' Computes the RMS of a numeric vector x.
#' @export
root_mean_square <- function(x) {
  squared <- x**2
  mean <- sum(squared)/length(squared)
  sqrt(mean)
}

#' Given a numeric vector x, returns a boolean vector whose ith truth value
#' indicates whether the ith element of x was above the 0.95 percentile of
#' rms(x).
#' @export
rms_spindle <- function(x, λ=0.95) {
  rms <- root_mean_square(x)
  threshold <- quantile(rms, λ)
  spindles <- sapply(rms, function(x) { x >= threshold })
  return(spindles)
}

#' @export
window_eeg_data <- function(data, window_size) {
    data$Windows <- cut(data$Time,
                             breaks = seq(0, max(data$Time), by = window_size),
                             include.lowest = TRUE)
    return(data)
}


#' Given a boolean vector x, sets to FALSE all sequences of TRUE values
#' with less than min or more than max elements.
#' In the context of spindle detection, the boolean vector is understood
#' signal which indexes of the EEG sample were found to contain spindles.
#' @export
remove_spindles_out_of_range <- function(x, min, max) {
  substitute <- c()
  run <- rle(x)
  for (i in 1:length(run$values)) {
    value <- run$values[i]
    length <- run$lengths[i]
    if (value == TRUE && length >= min && length <= max) {
        substitute <- c(substitute, rep(TRUE, length))
    }else{
      substitute <- c(substitute, rep(FALSE, length))
    }
  }
  return(substitute)
}

#' Sigma index method as described in O'Reilly & Nielsen (2015) and
#' Hupponen et al. (2007).
#' Parameters:
#' x: numeric vector
#' fs = 500: sampling frequency
#' 
#' Returns sigma index of vector x of EEG samples.
#' @export
sigma_index <- function(x, fs=500){
    amp_spectrum <- amplitude_spectrum(x, fs)
    spindle_band <- subset(amp_spectrum, Frequency >= 10.5 & Frequency <= 16.0)
    # mean_low has low alpha and theta bands.
    # mean_high has beta band.
    low <- subset(amp_spectrum, Frequency >= 4 & Frequency <= 10)
    high <- subset(amp_spectrum, Frequency >= 20 & Frequency <= 40)
    alpha <- subset(amp_spectrum, Frequency >= 7.5 & Frequency <= 10)
    mean_low <- mean(low$Amplitude)
    mean_high <- mean(high$Amplitude)
    max_alpha <- max(alpha$Amplitude)
    max <- max(spindle_band$Amplitude)
    if (max_alpha > max) {
        max <- 0
    }
    return(max / ((mean_low + mean_high) / 2))
}

#' @export
relative_spindle_power <- function(x, fs = 500) {
    amp_spectrum <- amplitude_spectrum(x, fs)
    spindle_band <- subset(amp_spectrum, Frequency >= 11 & Frequency <= 16.0)
    cross_band <- subset(amp_spectrum, Frequency >= 0.5 & Frequency <= 40)
    return(sum(spindle_band$Amplitude) / sum(cross_band$Amplitude))
}

#' @export
amplitude_spectrum <- function(x, fs = 500) {
    # Mean centering
    x <- x - mean(x)
    n <- length(x)
    # Apply Hanning window to the vector
    hanning_window <- gsignal::hann(n)
    x <- x * hanning_window
    #Zero-padd
    x <- c(x, rep(0, 512 - n))
    x <- fft(x)
    x <- (2 * abs(x)) / sum(hanning_window)
    Df <- fs / 512
    frequencies <- seq(0, (length(x) - 1)) * Df
    amp_spectrum <- data.frame(Frequency = frequencies, Amplitude = x)
    return(amp_spectrum)
}
