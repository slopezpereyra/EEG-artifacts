source("R/misc.R")

# Scripts on this file are dedicated to plotting EEG data
# and M-CAPA analysis results.

library(ggplot2)
library(lubridate)
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


#' @export
setGeneric(
    "plot_analysis_channel",
    function(object, chan, size) {
        standardGeneric("plot_analysis_channel")
    }
)


#' @export
setMethod(
    "plot_analysis_channel",
    "analysis",
    function(object, chan, size = 1) {
        canoms <- dplyr::filter(object@canoms, variate == chan)
        panoms <- dplyr::filter(object@panoms, variate == chan)
        data <- object@eeg@data

        # Get all times between start time and end time of each canom
        intervals <- mapply(function(x, y) data[x:y, 1], canoms$start, canoms$end)
        intervals <- union(intervals, panoms$Time)
        # Get all indexes between start and end of canoms
        locations <- mapply(function(x, y) x:y, canoms$start, canoms$end)
        # Unite with point anomalies
        locations <- union(unlist(locations), unlist(panoms$location))
        time_of_anomalies <- as_datetime(unlist(intervals))
        values <- unlist(data[locations, chan + 1])
        df <- tibble(A = time_of_anomalies, B = values)

        eeg <- plot_channel(object@eeg, channel = chan)
        p <- eeg +
            geom_point(
                data = df, aes(A, B),
                inherit.aes = FALSE, color = "red",
                size = size
            )
        return(p)
    }
)

#' @export
setGeneric(
    "plot_analysis",
    function(object, size) {
        standardGeneric("plot_analysis")
    }
)

#' Plots all channels with any kind of anomalies.
#'
#' @param object An Analysis object.
#' @return A plot_grid object.
#' @export
setMethod(
    "plot_analysis",
    "analysis",
    function(object, size = 1) {
        plots <- list()
        channels <- get_anomalous_channels(object)
        for (channel in channels)
        {
            p <- plot_analysis_channel(object, channel, size = size)
            plots[[channel]] <- p
        }
        return(plot_grid(plotlist = plots, align = "v", ncol = 1))
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


#' @export
setGeneric(
    "standardize_strengths",
    function(object, f) {
        standardGeneric("standardize_strengths")
    }
)


#' Normalizes point and collective anomalies' strengths
#' given a normalizing function f.
#'
#'
#' @return An analysis object.
#' @export
setMethod(
    "standardize_strengths",
    "analysis",
    function(object, f) {
        cstandard <- f(object@canoms$mean.change)
        pstandard <- f(object@panoms$strength)
        object@panoms$strength <- pstandard
        object@canoms$mean.change <- cstandard
        return(object)
    }
)
#' @export
setGeneric(
    "merge",
    function(object, an) {
        standardGeneric("merge")
    }
)


#' Merges two analysis into one
#'
#' @param object An analysis object
#' @param an The analysis object being merged into the first
#'
#' @return An analysis object.
#' @export
setMethod(
    "merge",
    "analysis",
    function(object, an) {
        canoms <- full_join(object@canoms, an@canoms)
        panoms <- full_join(object@panoms, an@panoms)

        data <- full_join(object@eeg@data, an@eeg@data)
        eeg <- new("eeg", data = data)

        return(new("analysis",
            canoms = canoms,
            panoms = panoms,
            eeg = eeg
        ))
    }
)

#' @export
setGeneric(
    "sfilter",
    function(object, x, f = minmax_normalization) {
        standardGeneric("sfilter")
    }
)


#' Filters an analysis object so as to keep only those
#' anomalies whose strength is greater than x.
#'
#' @param object An analysis object
#' @param x The strength threshold
#' @param f A function used to standardize strengths (otherwise collective and point
#' anomalies are incomparable). Defaults to min-max normalization.
#'
#' @return An analysis object.
#' @export
setMethod(
    "sfilter",
    "analysis",
    function(object, x, f = minmax_normalization) {
        standardized <- standardize_strengths(object, f)
        canoms <- standardized@canoms %>% dplyr::filter(mean.change >= x)
        panoms <- standardized@panoms %>% dplyr::filter(strength >= x)

        return(new("analysis",
            canoms = canoms,
            panoms = panoms,
            eeg = object@eeg
        ))
    }
)
