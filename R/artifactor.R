# Scripts on this file are dedicated to performing M-CAPA analysis
# on EEG data.

library(methods)
source("R/eeg.R")
source("R/cln.R")
source("R/analysis.r")


#' @export
setGeneric(
    "artf",
    function(eeg,
             s = -1, e = -1,
             res = 1,
             alpha = 8, beta = 1) {
        standardGeneric("artf")
    }
)

#' @export
setGeneric(
    "artf_stepwise",
    function(eeg,
             step_size = 30,
             res = 1,
             alpha = 8) {
        standardGeneric("artf_stepwise")
    }
)

#' Perform M-CAPA analysis on EEG data in a given range of seconds
#' and return results filtered by anomaly strength.
#'
#' @param eeg An eeg object
#' @param s int First second of timespan to analyze
#' @param e int Last second of timespan to analyze
#' @param res int Resolution at which to perform analysis
#' @param alpha float Threshold of strength significance for
#' collective anomalies
#' @param beta float Threshold of strength significance for point anomalies
#' @param time bool Show process time?
#' @return An analysis object
#' @export
#'
#' @export
setMethod(
    "artf",
    "eeg",
    function(eeg, s = -1, e = -1, res = 1, alpha = 8, beta = 1) {
        s <- ifelse(s != -1, s, eeg@data$Time[1])
        e <- ifelse(e != -1, e, tail(eeg@data$Time, n = 1))
        start_time <- Sys.time()
        eeg <- eeg %>%
            subset_eeg(s, e) %>%
            resample_eeg(res)
        analysis <- anomaly::capa.mv(eeg@data[-1], type = "mean")
        canoms <- anomaly::collective_anomalies(analysis) %>%
            dplyr::filter(mean.change >= alpha) %>%
            set_timevars(data = eeg@data) %>%
            tibble::as_tibble()
        panoms <- anomaly::point_anomalies(analysis) %>%
            dplyr::filter(strength >= beta) %>%
            set_timevars(data = eeg@data) %>%
            tibble::as_tibble()
        results <- new("analysis",
            canoms = canoms,
            panoms = panoms,
            eeg = eeg
        )
        end_time <- Sys.time()
        print(paste(
            "Analysis completed in ",
            (end_time - start_time), " seconds"
        ))
        return(results)
    }
)

#' Performs stepwise (or epoch by epoch) CAPA analysis on EEG data.
#'
#' @param eeg An eeg object.
#' @param step_size Size in seconds of each step.
#' @param alpha Threshold of strength significance for collective anomalies.
#' @return An analysis object.
#' @export
setMethod(
    "artf_stepwise",
    "eeg",
    function(eeg, step_size = 30, alpha = 8) {
        # Set epochs for grouping
        t <- set_epochs(eeg@data, epoch = step_size) %>% head(-1)
        mps <- get_sampling_frequency(eeg) * step_size # measures per step
        grouped <- dplyr::group_by(t[-1], Epoch) %>%
            dplyr::group_map(~ anomaly::capa.mv(x = .x, type = "mean"))
        canoms <- grouped %>% lapply(function(x) anomaly::collective_anomalies(x) %>% dplyr::filter(mean.change > 8))
        panoms <- grouped %>% lapply(function(x) anomaly::point_anomalies(x))

        canoms <- mapply(function(x, y) x %>% dplyr::mutate(start = start + mps * (y - 1), end = end + mps * (y - 1)),
            canoms, seq_along(canoms),
            SIMPLIFY = FALSE
        ) %>%
            dplyr::bind_rows() %>%
            set_timevars(eeg@data)

        panoms <- mapply(function(x, y) x %>% mutate(location = location + mps * (y - 1)),
            panoms, seq_along(panoms),
            SIMPLIFY = FALSE
        ) %>%
            dplyr::bind_rows() %>%
            set_timevars(eeg@data)

        an <- new("analysis", canoms = canoms, panoms = panoms, eeg = eeg)
        return(an)
    }
)

# Helper function
set_epochs <- function(df, epoch = 30) {
    # Get quotient and remainder of euclidean division
    # of each time in seconds by 30.
    q <- df$Time %/% epoch
    df <- df %>% tibble::add_column(Epoch = q + 1, .after = 1)
    return(df)
}
