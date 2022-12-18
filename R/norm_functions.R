
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
