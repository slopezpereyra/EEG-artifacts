
#' @export
setGeneric(
    "spectogram",
    function(object, channel, hcolors) standardGeneric("spectogram")
)

#' @export
setGeneric(
    "psd",
    function(object, epoch = 30) standardGeneric("psd")
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
    "psd_chan",
    "eeg",
    function(object, channel, epoch = 30) {
        t <- set_epochs(object@data, epoch) %>% head(-1)
        fs <- get_sampling_frequency(object)
        by_epoch <- tapply(
            unlist(t[-1][channel + 1]), # +1 accounts for epoch column
            t$Epoch,
            function(x) rsleep::pwelch(x, sRate = fs)
        )
        avg <- Reduce("+", lapply(by_epoch, "[[", "psd")) / length(by_epoch)
        df <- avg %>%
            tibble::as_tibble() %>%
            tibble::add_column(Fqc = by_epoch[[1]]$hz, .before = "value")
        return(df)
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
    function(object, epoch = 30) {
        t <- set_epochs(object@data, epoch) %>% head(-1)
        fs <- get_sampling_frequency(object)
        by_epoch <- dplyr::group_by(t[-1], Epoch) %>%
            # Cast pwelch PSD over each column on each each group.
            dplyr::group_map(~ apply(.x,
                FUN = function(y) rsleep::pwelch(y, sRate = fs),
                MARGIN = 2
            ))
        # Extract PSD of each channel over all epochs
        # (tibble 1 = epoch_1, ..., tibble_n = epoch_n)
        psds <- lapply(
            by_epoch,
            function(x) tibble::as_tibble(lapply(x, "[[", "psd"))
        )
        # Compute the average spectrum across epochs for each channel.
        avgs <- tibble::as_tibble(Reduce("+", psds) / length(by_epoch))
        # Add hz variable
        psd <- tibble::add_column(avgs,
            Fqc = by_epoch[[1]][[1]]$hz,
            .before = colnames(avgs)[1]
        )
        return(tibble::as_tibble(psd))
    }
)

#' Given an PSD data frame as returned by the psd(eeg, ...) function,
#' plot the spectrum of all channels.
#'
#' @param object An data frame as returned by the psd(eeg, ...) function.
#' @return A ggplot object.
#'
#' @export
plot_psd <- function(psd, xlim = 30) {
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
iplot_psd <- function(psd) {
    psd <- reshape2::melt(x, id.vars = "Fqc")
    fig <- plotly::plot_ly(
        psd,
        type = "scatter",
        mode = "lines"
    ) %>%
        plotly::add_trace(
            x = ~Fqc, y = ~value, color = ~variable
        )
    return(fig)
}
