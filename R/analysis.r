
# Scripts on this file are dedicated to plotting EEG data
# and M-CAPA analysis results.


#' Analysis class.
#' @slot canoms A collective anomalies data frame.
#' @slot signals A point anomalies data frame.
#' @slot The eeg upon which analysis is conducted.
#' @export
setClass("analysis",
    slots = list(
        canoms = "data.frame",
        panoms = "data.frame"
    )
)

# --------- Analysis Generics -------------

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
    function(object, chan, data) standardGeneric("set_plot_data")
)

#' @export
setGeneric(
    "plot_analysis_channel",
    function(object, chan, eeg, size = 0.2) standardGeneric("plot_analysis_channel")
)

#' @export
setGeneric(
    "plot_analysis",
    function(object, eeg, size = 0.2) standardGeneric("plot_analysis")
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
methods::setGeneric("extract_epochs", function(object) standardGeneric("extract_epochs"))


#' Show method for the EEG class that calls View on the canoms panoms
#' and eeg slots.
#' 
#' @param object An Analysis object.
#' @export
setMethod(
    "show",
    "analysis",
    function(object) {
        View(object@canoms)
        View(object@panoms)
    }
)

#' Check if either canoms or panoms are non-empty data frames.
#' @param object An analysis object.
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
#' anomalies in the analysis object.
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



#' Given an analysis object and an EEG channel, returns a data frame
#' whose format is idoneous for plotting the analysis results.
#' Such format consists of the EEG Time measure (x axis), the raw
#' EEG data of the chanel and a column that is filled with NaN except
#' at each time instance where any type of anomaly was found.
#' This method is not intended for isolated calls, but is implicitely used
#' in other plotting methods.
#'
#' @param object An Analysis object.
#' @param chan An integer that points to the EEG channel.
#' @param data Original EEG data for this analysis
#' @return A data frame
#' @export
setMethod(
    "set_plot_data",
    "analysis",
    function(object, chan, data) {
        canoms <- dplyr::filter(object@canoms, variate == chan)
        panoms <- dplyr::filter(object@panoms, variate == chan)
        # Get all indexes between start and end of canoms
        locations <- mapply(function(x, y) x:y, canoms$start, canoms$end)
        # Unite with point anomalies
        locations <- union(unlist(locations), unlist(panoms$location)) %>%
                    as.integer()
        time_of_anomalies <- lubridate::as_datetime(unlist(data[locations, 1]))
        values <- unlist(data[locations, chan + 1])
        df <- tibble::tibble(A = time_of_anomalies, B = values)
        return(df)
    }
)

#' Plots analysis results for a specific channel.
#'
#' @param object An Analysis object.
#' @param chan An integer that points to the EEG channel.
#' @param size (optional) The size of the red dots that signal an anomaly.
#' @return A data frame
#' @export
setMethod(
    "plot_analysis_channel",
    "analysis",
    function(object, chan, eeg, size = 0.2) {
        df <- set_plot_data(object, chan, eeg@data)
        eeg <- plot_channel(eeg, channel = chan)
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
#' @param x An Analysis object.
#' @return A plot_grid object.
#' @export
setMethod(
    "plot_analysis",
    "analysis",
    function(object, eeg, size = 0.2) {
        plots <- list()
        channels <- get_anomalous_channels(object)
        for (channel in channels) {
            p <- plot_analysis_channel(object, channel, eeg, size)
            plots[[channel]] <- p
        }
        return(cowplot::plot_grid(plotlist = plots, align = "v", ncol = 1))
    }
)


#' Normalizes the strengths of point and collective anomalies
#' given a normalizing function f. 
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


#' Filters an analysis object so as to keep only those
#' anomalies whose normalized strength is greater than x.
#'
#' @param object An analysis object
#' @param x The strength threshold
#' @param f A function used to normalize strengths (otherwise collective
#' and point anomalies are incomparable). Defaults to min-max normalization.
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
            panoms = panoms
        ))
    }
)

#' Given an analysis objects, returns a data frame with the epoch-subepoch pairs
#' found to contain anomalies, and their respective average strength.
#'
#' @param object An analysis object.
#' @return A data frame containing epoch-subepoch pairs
#' found to contain artifacts as well as the average
#' strength of the artifacts for each subepoch.
#' @export
methods::setMethod(
    "extract_epochs",
    "analysis",
    function(object) {
        a <- object@canoms %>%
            dplyr::group_by(Epoch, Subepoch) %>%
            dplyr::summarise_at(vars(mean.change), list(cAnomStrength = mean))
        b <- object@panoms %>%
            dplyr::group_by(Epoch, Subepoch) %>%
            dplyr::summarise_at(vars(strength), list(pAnomStrength = mean))
        return(dplyr::full_join(a, b))
    }
)
