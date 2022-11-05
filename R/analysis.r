

# Scripts on this file are dedicated to plotting EEG data
# and M-CAPA analysis results.

library(ggplot2)
library(cowplot)
source("R/eeg.R")

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
    "has_anoms",
    function(object, step, ...) {
        standardGeneric("has_anoms")
    }
)

#' Check if either canoms or panoms are non-empty data frames.
#' @param object An Analysis object.
#'
#' @return bool
#' @export
setMethod(
    "has_anoms",
    "analysis",
    function(object) {
        return(nrow(object@canoms) > 0 | nrow(object@panoms) > 0)
    }
)

#' @export
setGeneric(
    "get_anomalous_channels",
    function(object, channel) {
        standardGeneric("get_anomalous_channels")
    }
)


#' Getter function retrieving all EEG channels containing
#' anomalies in the Analysis object.
#'
#' @param object An Analysis object.
#' @return A vector containing all anomalous channels in the analysis.
#' @export
setMethod(
    "get_anomalous_channels",
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
    "plot_clusters",
    function(object, channel) {
        standardGeneric("plot_clusters")
    }
)

#' Plots only collective anomalies in one of the EEG channels.
#'
#' @param object An Analysis object.
#' @param channel An integer.
#' @return A ggplot.
setMethod(
    "plot_clusters",
    "analysis",
    function(object, channel) {
        canoms <- object@canoms %>%
            dplyr::filter(variate == channel)
        areas <- tibble(
            xmin = as_datetime(canoms$Time),
            xmax = as_datetime(object@eeg@data$Time[canoms$end]),
            ymin = -315, ymax = 315, alpha = 0.3
        )
        eeg <- plot_channel(object@eeg, channel = channel)
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
    "plot_points",
    function(object, channel) {
        standardGeneric("plot_points")
    }
)


#' Plots only point anomalies in one of the EEG channels.
#'
#' @param object An Analysis object.
#' @param channel An integer.
#' @return A ggplot.
#' @export
setMethod(
    "plot_points",
    "analysis",
    function(object, channel) {
        points <- object@panoms %>%
            dplyr::filter(variate == channel)
        point_values <- unlist(object@eeg@data[-1][points$location, channel])
        point_values <- unlist(object@eeg@data[-1][points$location, channel])
        points <- tibble(A = as_datetime(points$Time), B = point_values)

        eeg <- plot_channel(object@eeg, channel = channel)
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
    "plot_anoms",
    function(object, channel) {
        standardGeneric("plot_anoms")
    }
)


#' Plots both point and collective anomalies in one of the EEG channels.
#'
#' @param object An Analysis object.
#' @param channel An integer.
#' @return A ggplot.
#' @export
setMethod(
    "plot_anoms",
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

        eeg <- plot_channel(object@eeg, channel = channel)
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
    "plot_analyzed_channel",
    function(object, channel) {
        standardGeneric("plot_analyzed_channel")
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
    "plot_analyzed_channel",
    "analysis",
    function(object, channel) {
        canoms <- object@canoms %>% dplyr::filter(variate == channel)
        point_locs <- object@panoms %>%
            dplyr::filter(variate == channel) %>%
            .$location
        if (nrow(canoms) == 0 & length(point_locs) == 0) {
            return(plot.egg(object, channel))
        } else if (nrow(canoms) == 0) {
            return(plot_points(object, channel))
        } else {
            return(plot_anoms(object, channel))
        }
        return(plot)
    }
)

#' @export
setGeneric(
    "plot_analyzed_channels",
    function(object) {
        standardGeneric("plot_analyzed_channels")
    }
)

#' Plots all channels with any kind of anomalies.
#'
#' @param object An Analysis object.
#' @return A plot_grid object.
#' @export
setMethod(
    "plot_analyzed_channels",
    "analysis",
    function(object) {
        plots <- list()
        channels <- get_anomalous_channels(object)
        for (channel in channels)
        {
            p <- plot_analyzed_channel(object, channel)
            plots[[channel]] <- p
        }
        return(plot_grid(plotlist = plots, align = "v", ncol = 1))
    }
)

#' @export
setGeneric(
    "plot_analysis",
    function(object, save = FALSE) {
        standardGeneric("plot_analysis")
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
    "plot_analysis",
    "analysis",
    function(object, save = FALSE) {
        anom_plot <- plot_analyzed_channels(object)
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



#' @export
setGeneric(
    "set_chan_for_iplot",
    function(object, chan) {
        standardGeneric("set_chan_for_iplot")
    }
)

#' Formats a channel's analysis results to a data frame
#' specifically designed for interactive plotting using
#' Python's plotly library.
#'
#' @param object An Analysis object.
#' @param chan The channel to be formatted for intercative plotting.
#'
#' @return A tibble.
#' @export
setMethod(
    "set_chan_for_iplot",
    "analysis",
    function(object, chan) {
        df <- object@eeg@data[, c(1, (chan + 1))]
        df[, c("anoms", "strength")] <- NA
        canoms <- object@canoms %>% dplyr::filter(variate == chan)
        panoms <- object@panoms %>% dplyr::filter(variate == chan)
        # Set point anomalies appropriately
        point_locations <- unlist(panoms$location)
        point_values <- unlist(df[panoms$location, 2])
        point_strengths <- unlist(panoms$strength)
        df$anoms[point_locations] <- point_values
        df$strength[point_locations] <- point_strengths
        # Set channel anomalies
        canoms_positions <- list()
        canoms_strengths <- list()
        for (i in seq_len(nrow(canoms))) {
            if (canoms[i, ][["end"]] > nrow(df)) {
                break
            }
            start_end <- canoms[i, ][["start"]]:canoms[i, ][["end"]]
            strengths <- rep(canoms[i, ][["mean.change"]], length(start_end))
            canoms_positions <- append(canoms_positions, start_end)
            canoms_strengths <- append(canoms_strengths, strengths)
        }
        canoms_positions <- unlist(canoms_positions)
        canoms_strengths <- unlist(canoms_strengths)
        canom_values <- unlist(df[canoms_positions, 2])

        # canom_values <- df[canoms_positions, 2] %>%
        #    as.list() %>%
        df$anoms[canoms_positions] <- canom_values
        df$strength[canoms_positions] <- canoms_strengths
        return(df)
    }
)


#' @export
setGeneric(
    "set_for_iplot",
    function(object, save = FALSE) {
        standardGeneric("set_for_iplot")
    }
)


#' Sets the necessary data frame for interactive analysis plotting using
#' Python's plotly library.
#'
#' @param object An Analysis object.
#' @param save A bool declaring whether to save the resulting
#' data frame in in .csv format or not.
#'
#' @return A tibble.
#' @export
setMethod(
    "set_for_iplot",
    "analysis",
    function(object, save = FALSE) {
        nchans <- length(object@eeg@data[-1])
        df <- set_chan_for_iplot(object, 1)
        for (chan in 2:nchans) {
            iplot_data <- set_chan_for_iplot(object, chan)
            iplot_data$Time <- NULL
            df <- df %>% bind_cols(iplot_data)
        }
        if (save == TRUE) {
            write_csv(df, "iplot_data.csv")
        }
        return(df)
    }
)
