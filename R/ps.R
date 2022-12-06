source("R/artifactor.R")
library(reshape2)
library(rsleep)
library(gsignal)


#' @export
setGeneric("spectogram", function(object, channel, hcolors) standardGeneric("spectogram"))

#' @export
setGeneric("psd", function(object, method = "welch") standardGeneric("psd"))

#' @export
setGeneric("get_channel_psd", function(object, channel) standardGeneric("get_channel_psd"))


#' Produces the spectogram of an EEG object.
#' @param object The EEG object to be analyzed.
#' @param channel The channel to be analyzed.
setMethod(
    "spectogram",
    "eeg",
    function(object, channel, hcolors = 10) {
        fs <- get_sampling_frequency(object)
        x <- gsignal::specgram(unlist(object@data[channel + 1]), fs = fs)
        plot(x, col = grDevices::heat.colors(hcolors))
    }
)
#' Given an eeg object and a channel's column index,
#' compute the power spectrum density of  the channel
#' using Welch's method.
#'
#' @param object An eeg object.
#' @param channel A column number
#'
#' @return An data frame containing the power density
#' for every frequency composing the given channel's signal.
#'
#' @export
setMethod(
    "get_channel_psd",
    "eeg",
    function(object, channel) {
        vec <- unlist(object@data[channel])
        fs <- get_sampling_frequency(object)
        periodogram <- rsleep::pwelch(vec, sRate = fs)
        return(as_tibble(periodogram))
    }
)

#' Given an eeg object, compute the power spectrum of every
#' channel.
#'
#' @param object An eeg object.
#'
#' @return A data frame with the spectrum density for every frequency
#' for every channel in the  EEG record.
#'
#' @export
setMethod(
    "psd",
    "eeg",
    function(object) {
        df <- get_channel_psd(object, 2)
        for (chan in 3:(ncol(object@data) - 1)) {
            df <- add_column(df, get_channel_psd(object, chan)[2], .name_repair = "unique")
        }
        names <- colnames(object@data)
        names[1] <- "Fqc"
        colnames(df) <- names
        return(df)
    }
)

#' @export
plot_psd <- function(psd, xlim = 30) {
    tall_format <- melt(psd, id.vars = "Fqc")
    p <- ggplot(tall_format, aes(Fqc, value, col = variable)) +
        geom_line() +
        xlim(c(0, xlim))
    #    geom_vline(xintercept = 4, linetype = "dashed") +
    #    geom_vline(xintercept = 7, linetype = "dashed") +
    #    geom_vline(xintercept = 12, linetype = "dashed") +
    #    geom_vline(xintercept = 30, linetype = "dashed")
    p
}
