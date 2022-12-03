

znormalization <- function(x) {
    return((x - mean(x)) / sd(x))
}

get_interval <- function(x, y) {
    return(x:y)
}
