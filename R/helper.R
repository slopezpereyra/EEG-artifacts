
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
    x  <- x %>% dplyr::select(dplyr::starts_with("EEG") | dplyr::starts_with("EOG")) %>%
                tidyr::unnest(cols = dplyr::everything()) %>%
                dplyr::mutate(Time = unlist(t)) %>%
                dplyr::relocate(Time, .before = 1)

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

