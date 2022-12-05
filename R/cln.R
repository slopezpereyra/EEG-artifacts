

#' Given a canoms or panoms data frame, returns a vector containing
#' the time of occurrence of each anomaly on the EEG record.
#'
#' @param anoms A canoms or panoms data frame.
#'
#' @return A numeric vector
get_time <- function(anoms, data) {
    return(data$Time[unlist(anoms[1])]) # Remove seconds_to_period()
}

#' Given a canoms or panoms data frame, returns an identical
#' canoms or panoms with the inclusion of epoch and subepoch columns.
#'
#' @param anoms A canoms or panoms data frame.
#'
#' @return A data frame
set_anom_epoch <- function(anoms) {
    if (nrow(anoms) == 0) {
        return(anoms)
    }

    # Get quotient and remainder of euclidean division
    # of each time in seconds by 30.
    q <- anoms$Time %/% 30
    r <- anoms$Time %% 30
    # Get quotient of euclidean division of r by 5
    qprime <- r %/% 5

    anoms$Epoch <- q + 1
    anoms$Subepoch <- qprime + 1
    return(anoms)
}

#' Given a canoms or panoms data frame and its origin, returns an identical
#' canoms or panoms data frame with the inclusion of all time variables:
#' Time, Epoch and Subepoch.
#'
#' @param anoms A canoms or panoms data frame.
#' @param data The data where anomalies were detected.
#' @return A data frame
set_timevars <- function(anoms, data) {
    if (nrow(anoms) == 0) {
        return(anoms)
    }
    anoms$Time <- get_time(anoms, data)
    anoms <- set_anom_epoch(anoms)
    return(anoms)
}
