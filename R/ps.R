
#' @export
setGeneric(
    "spectogram",
    function(object, channel, maxFreq = 30, freq = 4) standardGeneric("spectogram")
)

#' @export
setGeneric(
    "psd",
    function(object) standardGeneric("psd")
)
#' @export
setGeneric(
    "psd_chan",
    function(object, channel, epoch = 30) standardGeneric("psd_chan")
)


#' Produces the spectogram of an EEG object.
#' @param object The EEG object to be analyzed.
#' @param channel The channel to be analyzed.
setMethod(
    "spectogram",
    "eeg",
    function(object, channel, maxFreq = 30, freq = 4) {
        fs <- get_sampling_frequency(object)
        rsleep::spectrogram(unlist(object@data[-1][channel]), sRate = fs, maxFreq = maxFreq, freq = freq)
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
    "psd_chan",
    "eeg",
    function(object, channel) {
        fs <- get_sampling_frequency(object)
        pwelch <- gsignal::pwelch(as.matrix(object@data[-1][channel]), fs = fs)
        psd <- pwelch$spec %>%
            apply(log10, MARGIN = 2) %>%
            tibble::as_tibble()
        psd$Fqc <- pwelch$freq
        return(psd)
    }
)

#' Given an eeg object, compute the average power spectrum density
#' of each channel across epochs.
#'
#' @param object An eeg object.
#' @return A data frame with the spectrum density for every frequency
#' for every channel in the  EEG record.
#'
#' @export
setMethod(
    "psd",
    "eeg",
    function(object) {
        fs <- get_sampling_frequency(object)
        pwelch <- gsignal::pwelch(as.matrix(object@data[-1]), fs = fs)
        psd <- pwelch$spec %>%
            apply(log10, MARGIN = 2) %>%
            tibble::as_tibble()
        psd$Fqc <- pwelch$freq
        return(psd)
    }
)

#' Given an PSD data frame as returned by the psd(eeg, ...) function,
#' plot the spectrum of all channels.
#'
#' @param object An data frame as returned by the psd(eeg, ...) function.
#' @return A ggplot object.
#'
#' @export
plot_psd <- function(psd, xlim = 250) {
    tall_format <- reshape2::melt(psd, id.vars = "Fqc")
    p <- ggplot2::ggplot(tall_format, ggplot2::aes(Fqc, value, col = variable)) +
        ggplot2::geom_line() +
        ggplot2::xlim(c(0, xlim))
    #    geom_vline(xintercept = 4, linetype = "dashed") +
    #    geom_vline(xintercept = 7, linetype = "dashed") +
    #    geom_vline(xintercept = 12, linetype = "dashed") +
    #    geom_vline(xintercept = 30, linetype = "dashed")
    return(p)
}


#' Given an PSD data frame as returned by the psd(eeg, ...) function,
#' create an interactive plot of the spectrum of all channels.
#'
#' @param object An data frame as returned by the psd(eeg, ...) function.
#' @return A plotly figure.
#'
#' @export
iplot_psd <- function(psd, xlim = 250) {
    psd <- reshape2::melt(psd, id.vars = "Fqc")
    fig <- plotly::plot_ly(
        psd,
        type = "scatter",
        mode = "lines"
    ) %>%
        plotly::add_trace(
            x = ~Fqc, y = ~value, color = ~variable
        ) %>%
        layout(
            xaxis = list(
                title = "Frequency in Hz",
                zeroline = F,
                range = c(0, xlim)
            ),
            yaxis = list(
                title = "log10 Spectrum",
                zeroline = F
            )
        )
    return(fig)
}
