

# Scripts on this file are dedicated to plotting EEG data
# and M-CAPA analysis results.

library(ggplot2)
library(cowplot)

#' Analysis class.
#' @slot canoms A collective anomalies data frame.
#' @slot signals A point anomalies data frame.
#' @slot The eeg upon which analysis is conducted.
#' @export
setClass("analysis",
    slots = list(
        canoms = "data.frame",
        panoms = "data.frame",
        eeg = "eeg"
    )
)

#' Show method for the EEG class that calls View on the canoms panoms
#' and eeg slots.
#' @param object An Analysis object.
#' @export
setMethod(
    "show",
    "analysis",
    function(object) {
        View(object@canoms)
        View(object@panoms)
        View(object@eeg@data)
    }
)

#' @export
setGeneric(
    "has.anomalies",
    function(object, step, ...) {
        standardGeneric("has.anomalies")
    }
)

#' Check if either canoms or panoms are non-empty data frames.
#' @param object An Analysis object.
#'
#' @return bool
#' @export
setMethod(
    "has.anomalies",
    "analysis",
    function(object) {
        return(nrow(object@canoms) > 0 | nrow(object@panoms) > 0)
    }
)

#' @export
setGeneric(
    "get.anomalous.channels",
    function(object, channel) {
        standardGeneric("get.anomalous.channels")
    }
)


#' Getter function retrieving all EEG channels containing
#' anomalies in the Analysis object.
#'
#' @param object An Analysis object.
#' @return A vector containing all anomalous channels in the analysis.
#' @export
setMethod(
    "get.anomalous.channels",
    "analysis",
    function(object) {
        chans <- union(
            unique(object@canoms$variate),
            unique(object@panoms$variate)
        )
        return(chans)
    }
)

setGeneric(
    "plot.clusters",
    function(object, channel) {
        standardGeneric("plot.clusters")
    }
)

#' Plots only collective anomalies in one of the EEG channels.
#'
#' @param object An Analysis object.
#' @param channel An integer.
#' @return A ggplot.
setMethod(
    "plot.clusters",
    "analysis",
    function(object, channel) {
        canoms <- object@canoms %>%
            dplyr::filter(variate == channel)
        areas <- tibble(
            xmin = as_datetime(canoms$Time),
            xmax = as_datetime(object@eeg@data$Time[canoms$end]),
            ymin = -315, ymax = 315, alpha = 0.3
        )
        eeg <- plot.channel(object@eeg, channel = channel)
        plot <- eeg +
            geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = 0.3, fill = "red",
                data = transform(areas, as.character(seq_len(nrow(areas)))),
                inherit.aes = FALSE
            ) + xlab("") + ylab(colnames(object@eeg@data[-1][channel]))
        return(plot)
    }
)

#' @export
setGeneric(
    "plot.points",
    function(object, channel) {
        standardGeneric("plot.points")
    }
)


#' Plots only point anomalies in one of the EEG channels.
#'
#' @param object An Analysis object.
#' @param channel An integer.
#' @return A ggplot.
#' @export
setMethod(
    "plot.points",
    "analysis",
    function(object, channel) {
        points <- object@panoms %>%
            dplyr::filter(variate == channel)
        point_values <- unlist(object@eeg@data[-1][points$location, channel])
        points <- tibble(A = as_datetime(points$Time), B = point_values)

        eeg <- plot.channel(object@eeg, channel = channel)
        plot <- eeg +
            geom_point(
                data = points, aes(A, B),
                inherit.aes = FALSE, color = "red"
            ) +
            xlab("") + ylab(colnames(object@eeg@data[-1][channel]))

        return(plot)
    }
)

#' @export
setGeneric(
    "plot.anomalies",
    function(object, channel) {
        standardGeneric("plot.anomalies")
    }
)


#' Plots both point and collective anomalies in one of the EEG channels.
#'
#' @param object An Analysis object.
#' @param channel An integer.
#' @return A ggplot.
#' @export
setMethod(
    "plot.anomalies",
    "analysis",
    function(object, channel) {
        points <- object@panoms %>%
            dplyr::filter(variate == channel)
        point_values <- unlist(object@eeg@data[-1][points$location, channel])
        points <- tibble(A = as_datetime(points$Time), B = point_values)

        canoms <- object@canoms %>% dplyr::filter(variate == channel)
        areas <- tibble(
            xmin = as_datetime(canoms$Time),
            xmax = as_datetime(object@eeg@data$Time[canoms$end]),
            ymin = -300, ymax = 300, alpha = 0.3
        )

        eeg <- plot.channel(object@eeg, channel = channel)
        plot <- eeg +
            geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = 0.3, fill = "red",
                data = transform(areas, as.character(1:nrow(areas))),
                inherit.aes = FALSE
            ) +
            geom_point(
                data = points, aes(A, B),
                inherit.aes = FALSE, color = "red"
            ) +
            xlab("") + ylab(colnames(object@eeg@data[-1][channel]))

        return(plot)
    }
)


#' @export
setGeneric(
    "plot.channel.analysis",
    function(object, channel) {
        standardGeneric("plot.channel.analysis")
    }
)

#' Wrapper function that plots an analyzed channel with
#' any kind of anomalies in it.
#'
#' @param object An Analysis object.
#' @param channel An integer.
#' @return A ggplot.
#' @export
setMethod(
    "plot.channel.analysis",
    "analysis",
    function(object, channel) {
        canoms <- object@canoms %>% dplyr::filter(variate == channel)
        point_locs <- object@panoms %>%
            dplyr::filter(variate == channel) %>%
            .$location
        if (nrow(canoms) == 0 & length(point_locs) == 0) {
            return(plot.egg(object, channel))
        } else if (nrow(canoms) == 0) {
            return(plot.points(object, channel))
        } else {
            return(plot.anomalies(object, channel))
        }
        return(plot)
    }
)

#' @export
setGeneric(
    "plot.channels.analysis",
    function(object) {
        standardGeneric("plot.channels.analysis")
    }
)

#' Plots all channels with any kind of anomalies.
#'
#' @param object An Analysis object.
#' @return A plot_grid object.
#' @export
setMethod(
    "plot.channels.analysis",
    "analysis",
    function(object) {
        plots <- list()
        channels <- get.anomalous.channels(object)
        for (channel in channels)
        {
            p <- plot.channel.analysis(object, channel)
            plots[[channel]] <- p
        }
        return(plot_grid(plotlist = plots, align = "v", ncol = 1))
    }
)

#' @export
setGeneric(
    "plot.analysis",
    function(object, save = FALSE) {
        standardGeneric("plot.analysis")
    }
)

#' Wrapper and standard analysis plotting function.
#' Plots all analyzed channels with the option of saving
#' the plot.
#'
#' @param object An Analysis object.
#' @param save A bool declaring whether to save the image in .png format or not.
#'
#' @return A plot_grid object.
#' @export
setMethod(
    "plot.analysis",
    "analysis",
    function(object, save = FALSE) {
        anom_plot <- plot.channels.analysis(object)
        if (save == TRUE) {
            dir.create("results")
            start <- seconds_to_period(object@eeg@data$Time[1]) %>% format.time()
            end <- seconds_to_period(tail(object@eeg@data$Time, 1)) %>% format.time()
            time <- paste(start, end, sep = " to ")
            ggsave(paste("results/", time, ".png", sep = ""), plot = anom_plot)
        }
        return(anom_plot)
    }
)
