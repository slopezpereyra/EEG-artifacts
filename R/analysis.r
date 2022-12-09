source("R/norm_functions.R")

# Scripts on this file are dedicated to plotting EEG data
# and M-CAPA analysis results.

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

# --------- Analysis Generics -------------

#' @export
setGeneric(
    "set_chan_for_iplot",
    function(object, chan) standardGeneric("set_chan_for_iplot")
)

#' @export
setGeneric(
    "has_anoms",
    function(object) standardGeneric("has_anoms")
)

#' @export
setGeneric(
    "get_anomalous_channels",
    function(object, channel) standardGeneric("get_anomalous_channels")
)

#' @export
setGeneric(
    "set_plot_data",
    function(object, chan) standardGeneric("set_plot_data")
)

#' @export
setGeneric(
    "plot_analysis_channel",
    function(object, chan, size = 0.2) standardGeneric("plot_analysis_channel")
)

#' @export
setGeneric(
    "plot_analysis",
    function(object, size = 0.2) standardGeneric("plot_analysis")
)

#' @export
setGeneric(
    "set_for_iplot",
    function(object, save = FALSE) standardGeneric("set_for_iplot")
)

#' @export
setGeneric(
    "standardize_strengths",
    function(object, f) standardGeneric("standardize_strengths")
)

#' @export
setGeneric(
    "sfilter",
    function(object, x, f = minmax_normalization) standardGeneric("sfilter")
)

#' @export
setGeneric(
    "merge",
    function(object, an) standardGeneric("merge")
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
setMethod(
    "set_plot_data",
    "analysis",
    function(object, chan) {
        canoms <- dplyr::filter(object@canoms, variate == chan)
        panoms <- dplyr::filter(object@panoms, variate == chan)
        data <- object@eeg@data
        # Get all indexes between start and end of canoms
        locations <- mapply(function(x, y) x:y, canoms$start, canoms$end)
        # Unite with point anomalies
        locations <- union(unlist(locations), unlist(panoms$location))
        time_of_anomalies <- lubridate::as_datetime(unlist(data[locations, 1]))
        values <- unlist(data[locations, chan + 1])
        df <- tibble::tibble(A = time_of_anomalies, B = values)
        return(df)
    }
)

#' @export
setMethod(
    "plot_analysis_channel",
    "analysis",
    function(object, chan, size = 0.2) {
        df <- set_plot_data(object, chan)

        eeg <- plot_channel(object@eeg, channel = chan)
        p <- eeg +
            ggplot2::geom_point(
                data = df, ggplot2::aes(A, B),
                inherit.aes = FALSE, color = "red",
                size = size
            )
        return(p)
    }
)

#' Plots all channels with any kind of anomalies.
#'
#' @param object An Analysis object.
#' @return A plot_grid object.
#' @export
setMethod(
    "plot",
    "analysis",
    function(x) {
        plots <- list()
        channels <- get_anomalous_channels(x)
        for (channel in channels) {
            p <- plot_analysis_channel(x, channel, size = 0.2)
            plots[[channel]] <- p
        }
        return(cowplot::plot_grid(plotlist = plots, align = "v", ncol = 1))
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
        df$anoms[canoms_positions] <- canom_values
        df$strength[canoms_positions] <- canoms_strengths
        return(df)
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
            df <- df %>% dplyr::bind_cols(iplot_data)
        }
        if (save == TRUE) {
            readr::write_csv(df, "iplot_data.csv")
        }
        return(df)
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
        canoms <- dplyr::full_join(object@canoms, an@canoms)
        panoms <- dplyr::full_join(object@panoms, an@panoms)

        data <- dplyr::full_join(object@eeg@data, an@eeg@data)
        eeg <- new("eeg", data = data)

        return(new("analysis",
            canoms = canoms,
            panoms = panoms,
            eeg = eeg
        ))
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
