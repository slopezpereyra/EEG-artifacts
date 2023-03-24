
# Scripts on this file are dedicated to plotting EEG data
# and M-CAPA analysis results.


# --------- Analysis Generics -------------

#' @export
setGeneric(
    "has_artifacts",
    function(eeg) standardGeneric("has_artifacts")
)

#' @export
setGeneric(
    "get_contaminated_channels",
    function(eeg, channel) standardGeneric("get_contaminated_channels")
)

#' @export
setGeneric(
    "set_plot_data",
    function(eeg, chan, data) standardGeneric("set_plot_data")
)

#' @export
setGeneric(
    "plot_channel_artifacts",
    function(eeg, chan, size = 0.2) standardGeneric("plot_channel_artifacts")
)

#' @export
setGeneric(
    "plot_artifacts",
    function(eeg, size = 0.2) standardGeneric("plot_artifacts")
)


#' @export
setGeneric(
    "standardize_strengths",
    function(eeg, f) standardGeneric("standardize_strengths")
)

#' @export
setGeneric(
    "sfilter",
    function(eeg, x, f = minmax_normalization) standardGeneric("sfilter")
)


#' @export
methods::setGeneric("extract_epochs", function(eeg) standardGeneric("extract_epochs"))


#' Check if either canoms or panoms are non-empty data frames.
#' @param eeg An analysis object.
#'
#' @return bool
#' @export
setMethod(
    "has_artifacts",
    "eeg",
    function(eeg) {
        return(nrow(eeg@canoms) > 0 | nrow(eeg@panoms) > 0)
    }
)


#' Getter function retrieving all EEG channels containing
#' anomalies in the analysis eeg.
#'
#' @param eeg An Analysis object.
#' @return A vector containing all anomalous channels in the analysis.
#' @export
setMethod(
    "get_contaminated_channels",
    "eeg",
    function(eeg) {
        chans <- union(
            unique(eeg@canoms$variate),
            unique(eeg@panoms$variate)
        )
        return(chans)
    }
)



#' Given an analysis eeg and an EEG channel, returns a data frame
#' whose format is idoneous for plotting the analysis results.
#' Such format consists of the EEG Time measure (x axis), the raw
#' EEG data of the chanel and a column that is filled with NaN except
#' at each time instance where any type of anomaly was found.
#' This method is not intended for isolated calls, but is implicitely used
#' in other plotting methods.
#'
#' @param eeg An Analysis object.
#' @param chan An integer that points to the EEG channel.
#' @param data Original EEG data for this analysis
#' @return A data frame
#' @export
setMethod(
    "set_plot_data",
    "eeg",
    function(eeg, chan, data) {
        canoms <- dplyr::filter(eeg@canoms, variate == chan)
        panoms <- dplyr::filter(eeg@panoms, variate == chan)
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
#' @param eeg An Analysis object.
#' @param chan An integer that points to the EEG channel.
#' @param size (optional) The size of the red dots that signal an anomaly.
#' @return A data frame
#' @export
setMethod(
    "plot_channel_artifacts",
    "eeg",
    function(eeg, chan, size = 0.2) {
        df <- set_plot_data(eeg, chan, object@data)
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
#' @param x An Analysis eeg.
#' @return A plot_grid eeg.
#' @export
setMethod(
    "plot_artifacts",
    "eeg",
    function(eeg, size = 0.2) {
        plots <- list()
        channels <- get_contaminated_channels(eeg)
        for (channel in channels) {
            p <- plot_channel_artifacts(eeg, channel, size)
            plots[[channel]] <- p
        }
        return(cowplot::plot_grid(plotlist = plots, align = "v", ncol = 1))
    }
)


#' Normalizes the strengths of point and collective anomalies
#' given a normalizing function f. 
#'
#' @return An analysis eeg.
#' @export
setMethod(
    "standardize_strengths",
    "eeg",
    function(eeg, f) {
        cstandard <- f(eeg@canoms$mean.change)
        pstandard <- f(eeg@panoms$strength)
        eeg@panoms$strength <- pstandard
        eeg@canoms$mean.change <- cstandard
        return(eeg)
    }
)


#' Filters an analysis eeg so as to keep only those
#' anomalies whose normalized strength is greater than x.
#'
#' @param eeg An analysis object
#' @param x The strength threshold
#' @param f A function used to normalize strengths (otherwise collective
#' and point anomalies are incomparable). Defaults to min-max normalization.
#'
#' @return An analysis eeg.
#' @export
setMethod(
    "sfilter",
    "eeg",
    function(eeg, x, f = minmax_normalization) {
        standardized <- standardize_strengths(eeg, f)
        canoms <- standardized@canoms %>% dplyr::filter(mean.change >= x)
        panoms <- standardized@panoms %>% dplyr::filter(strength >= x)

        return(new("eeg",
            canoms = canoms,
            panoms = panoms
        ))
    }
)

#' Helper function
#' @export
canoms_avg_epoch_strength <- function(df) {
    a <- df %>%
            dplyr::group_by(Epoch, Subepoch) %>%
            dplyr::summarise_at(vars(mean.change), list(cAnomStrength = mean))
    print(a)
    return(a)
}

#' Helper function
#' @export
panoms_avg_epoch_strength <- function(df) {
    a <- df %>%
            dplyr::group_by(Epoch, Subepoch) %>%
            dplyr::summarise_at(vars(strength), list(pAnomStrength = mean))
    return(a)

}

#' Given an analysis eeg, returns a data frame with the epoch-subepoch pairs
#' found to contain anomalies, and their respective average strength.
#'
#' @param eeg An analysis object.
#' @return A data frame containing epoch-subepoch pairs
#' found to contain artifacts as well as the average
#' strength of the artifacts for each subepoch.
#' @export
methods::setMethod(
    "extract_epochs",
    "eeg",
    function(eeg) {
        l <- list(eeg@canoms, eeg@panoms)
        i <- which(lapply(l, nrow) == 0)
        # If both data frames have length 0
        if (length(i) > 1) {
            stop("Empty analysis error")
        }
        # If both data frames are non-empty
        if (identical(i, integer(0))) {
            a <- canoms_avg_epoch_strength(eeg@canoms)
            b <- panoms_avg_epoch_strength(eeg@panoms)
            return(dplyr::full_join(a, b))
        }
        # If only one is non-empty
        if (i == 1) {
            return(panoms_avg_epoch_strength(eeg@panoms))
        }
        return(canoms_avg_epoch_strength(eeg@canoms))
    }
)
