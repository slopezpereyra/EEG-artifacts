
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
