#' Scripts on this file are devoted to the EEG class and its functions.

#' EEG class.
#' @slot data A data frame containing EEG records
#' @slot signals A data frame containing signal information (may be empty).
#' @export
setClass("eeg",
    slots = list(
        data = "data.frame",
        signals = "data.frame",
        canoms = "data.frame",
        panoms = 'data.frame'
    )
)

#' Show method for the EEG class that prints the data slot.
#' @param object An EEG object.
#' @export
methods::setMethod(
    "show",
    "eeg",
    function(object) {
        print(object@data)
    }
)

# ---------- EEG Generics ------------

#' @export
methods::setGeneric(
    "plot_channel",
    function(object, channel) standardGeneric("plot_channel")
)

#' @export
methods::setGeneric(
    "iplot",
    function(eeg) standardGeneric("iplot")
)

#' @export
methods::setGeneric(
    "bandpass",
    function(object, l, h) standardGeneric("bandpass")
)

#' @export
methods::setGeneric(
    "high_pass",
    function(object, n) standardGeneric("high_pass")
)

#' @export
methods::setGeneric(
    "subset_eeg",
    function(object, s, e) standardGeneric("subset_eeg")
)

#' @export
methods::setGeneric(
    "resample_eeg",
    function(object, n) standardGeneric("resample_eeg")
)

#' @export
methods::setGeneric(
    "get_epoch",
    function(object, which_epoch, epoch = 30) standardGeneric("get_epoch")
)

#' @export
methods::setGeneric(
    "get_sampling_frequency",
    function(object) standardGeneric("get_sampling_frequency")
)

#' @export
methods::setGeneric(
    "low_pass",
    function(object, n) standardGeneric("low_pass")
)

#' @export
methods::setGeneric(
    "remove_epoch",
    function(object, epoch) standardGeneric("remove_epoch")
)

#' @export
methods::setGeneric(
    "drop_epochs",
    function(object, epochs, epoch = 30) standardGeneric("drop_epochs")
)

#' @export
methods::setGeneric(
    "drop_subepochs",
    function(object, epochs, subepochs, epoch = 30) standardGeneric("drop_subepochs")
)

#' @export
methods::setGeneric(
    "artf_reject",
    function(object, analysis) standardGeneric("artf_reject")
)

#' Reads a .csv data file containing EEG data and an optional
#' signals file and return an eeg object. If the signals file is provided,
#' channel names are appropriately set.
#' If no signals file is provided, the method call is equivalent to
#' to simply using readr::read_csv(data_file)
#'
#' @param data_file .csv file containing eeg data.
#' @param signals_file .csv file containing signal information
#'
#' @return An eeg object.
#' @export
load_eeg <- function(data_file, signals_file = NULL) {
    data <- readr::read_csv(data_file)
    if (!is.null(signals_file)) {
        signals <- readr::read_csv(signals_file)
        colnames(data)[-1] <- signals$Label %>%
            stringr::str_remove("EEG ") %>%
            stringr::str_remove("EOG")
    } else {
        signals <- tibble::tibble()
    }

    return(new("eeg", data = data, signals = signals))
}


#' Subsets an EEG object from second s to e and
#' returns new EEG object whose data is that subset.
#'
#' @param object An eeg object.
#' @param s Starting time of the subset in seconds.
#' @param e Ending time of the subset in seconds.
#'
#' @return A new eeg whose data is the subset ranging from
#' second s to e of the object's data attribute.
#' @export
methods::setMethod(
    "subset_eeg",
    "eeg",
    function(object, s, e) {
        s_ind <- which(object@data$Time == s)
        e_ind <- which(object@data$Time == e)

        return(new("eeg",
            data = object@data[s_ind:e_ind, ],
            signals = object@signals
        ))
    }
)



#' Return a new EEG object containing only one every
#' n samples from the input EEG.
#'
#' @param object An eeg object.
#' @param n An integer.
#'
#' @return An EEG object
#' @export
methods::setMethod(
    "resample_eeg",
    "eeg",
    function(object, n) {
        df <- object@data[seq(1, nrow(object@data), n), ]
        return(new("eeg", data = df, signals = object@signals))
    }
)

#' Returns a new EEG object made up of only the specified epoch
#' of the input EEG.
#'
#' This function wraps a specific call of subset_eeg.
#'
#' @param object An eeg object.
#' @param which_epoch What epoch to subset.
#' @param epoch How many seconds make up an epoch? Defaults to 30.
#'
#' @return A new eeg whose data is the requested epoch.
#' @export
methods::setMethod(
    "get_epoch",
    "eeg",
    function(object, which_epoch, epoch = 30) {
        subset_eeg(object, epoch * (which_epoch - 1), epoch * (which_epoch))
    }
)


#' Returns the sampling frequency of the EEG object.
#'
#' @param object An eeg object.
#'
#' @return A number, the sampling frequency
#' @export
methods::setMethod(
    "get_sampling_frequency",
    "eeg",
    function(object) {
        delta_t <- object@data$Time[2] - object@data$Time[1]
        return(1 / delta_t)
    }
)

#' Wrapper function for low-pass filtering
#' a vector given a filtering frequency and
#' a sampling frequency.
#'
#' @export
vlow_pass <- function(vec, n, fs) {
    wpass <- n / (fs / 2) # Nyquist
    but <- gsignal::butter(5, wpass, "low", output = "Sos")
    low_pass <- gsignal::filter(but, unlist(vec))
    return(low_pass)
}

#' Wrapper function for high-pass filtering
#' a vector given a filtering frequency and
#' a sampling frequency.
#'
#' @export
vhigh_pass <- function(vec, n, fs) {
    wpass <- n / (fs / 2) # Nyquist
    but <- gsignal::butter(5, wpass, "high", output = "Sos")
    high_pass <- gsignal::filter(but, unlist(vec))
    return(high_pass)
}

#' Wrapper function for bandpass filtering
#' a vector given filtering frequencies and
#' a sampling frequency.
#'
#' @export
vbandpass <- function(vec, l, h, fs) {
    fpass <- c(l, h)
    wpass <- fpass / (fs / 2) # Nyquist
    but <- gsignal::butter(5, wpass, "pass", output = "Sos")
    pass <- gsignal::filter(but, unlist(vec))
    return(pass)
}

#' Given an EEG object and a numeric frequency n,
#' applies a low-pass Butterworth filter to all channels.
#'
#' @param object An eeg object.
#' @param n Filter frequency in Hz
#'
#' @return A new filtered EEG object
#' @export
methods::setMethod(
    "low_pass",
    "eeg",
    function(object, n) {
        fs <- get_sampling_frequency(object)
        df <- object@data
        filt_df <- apply(df[-1],
            MARGIN = 2,
            FUN = function(x) vlow_pass(x, n, fs),
            simplify = FALSE
        ) %>%
            tibble::as_tibble() %>%
            tibble::add_column(Time = df$Time, .before = colnames(df)[2])
        return(new("eeg", data = filt_df, signals = object@signals))
    }
)


#' Given an EEG object and a numeric frequency n,
#' applies a low-pass Butterworth filter to all channels.
#'
#' @param object An eeg object.
#' @param n Filter frequency in Hz
#'
#' @return A new filtered EEG object
#' @export
methods::setMethod(
    "high_pass",
    "eeg",
    function(object, n) {
        fs <- get_sampling_frequency(object)
        df <- object@data
        filt_df <- apply(df[-1],
            MARGIN = 2,
            FUN = function(x) vhigh_pass(x, n, fs),
            simplify = FALSE
        ) %>%
            tibble::as_tibble() %>%
            tibble::add_column(Time = df$Time, .before = colnames(df)[2])
        return(new("eeg", data = filt_df, signals = object@signals))
    }
)


#' Given an eeg object and a numeric frequency n,
#' applies a bandpass filter to all channels.
#'
#' @param object An eeg object.
#' @param l lower bound in Hz
#' @param h higher bound in Hz
#'
#' @return A new filtered EEG object
#' @export
methods::setMethod(
    "bandpass",
    "eeg",
    function(object, l, h) {
        fs <- get_sampling_frequency(object)
        df <- object@data
        filt_df <- apply(df[-1],
            MARGIN = 2,
            FUN = function(x) vbandpass(x, l, h, fs),
            simplify = FALSE
        ) %>%
            tibble::as_tibble() %>%
            tibble::add_column(Time = df$Time, .before = colnames(df)[2])
        return(new("eeg", data = filt_df, signals = object@signals))
    }
)


#' Given an EEG object and a channel integer n,
#' plots EEG record of the nth channel.
#'
#' @param object An eeg object.
#' @param channel An integer indicating index of channel to plot.
#'
#' @return A ggplot object.
#' @export
methods::setMethod(
    "plot_channel",
    "eeg",
    function(object, channel) {
        y <- object@data[-1][channel]
        p <- object@data %>%
            dplyr::mutate(Time = lubridate::as_datetime(Time)) %>%
            ggplot2::ggplot(
                ggplot2::aes(
                    Time,
                    unlist(y)
                )
            ) +
            ggplot2::geom_line() +
            ggplot2::scale_x_datetime(date_labels = "%H:%M:%S") +
            ggplot2::xlab("") +
            ggplot2::ylab(colnames(object@data[-1][channel]))

        return(p)
    }
)

#' Given an EEG object, plot all the EEG channels
#' in a vertical single-column layout.
#'
#' @param object An eeg object.
#'
#' @return A plot_grid object.
#' @export
methods::setMethod(
    "plot",
    "eeg",
    function(x) {
        plots <- list()
        for (channel in 1:(ncol(x@data) - 1)) {
            p <- plot_channel(x, channel)
            plots[[channel]] <- p
        }

        return(cowplot::plot_grid(plotlist = plots, align = "v", ncol = 1))
    }
)

#' Given an eeg object and an epoch, returns a new EEG object with said
#' epoch removed. For practical reasons, the extremes of the epoch are kept.
#'
#' @param object An eeg object.
#' @param epoch An natural number
#' @return An EEG object.
#' @export
methods::setMethod(
    "remove_epoch",
    "eeg",
    function(object, epoch) {
        s_ind <- which(object@data$Time == 30 * (epoch - 1)) + 1
        e_ind <- which(object@data$Time == 30 * epoch) - 1
        df <- object@data[-c(s_ind:e_ind), ]
        return(new("eeg", data = df, signals = object@signals))
    }
)

#' Given an EEG, return its interactive visualization. This method should not
#' be called on very large EEG objects (e.g., above an hour) without
#' previous resampling. Interactive visualizations are computationally
#' expensive.
#'
#' @param object An eeg object.
#' @return A plotly figure.
#' @export
methods::setMethod(
    "iplot",
    "eeg",
    function(eeg) {
        plots <- list()
        plots[[1]] <- plotly::plot_ly(
            eeg@data,
            type = "scatter",
            mode = "lines"
        ) %>%
            plotly::add_trace(
                x = ~Time, y = unlist(eeg@data[, 2]),
                name = colnames(eeg@data)[2]
            )
        for (i in 3:length(eeg@data)) {
            fig <- plotly::plot_ly(
                eeg@data,
                type = "scatter",
                mode = "lines"
            ) %>%
                plotly::add_trace(
                    x = ~Time, y = unlist(eeg@data[, i]),
                    name = colnames(eeg@data)[i]
                )
            plots[[i - 1]] <- fig
        }
        return(plotly::subplot(
            plots,
            nrows = length(plots),
            shareX = TRUE
        ))
    }
)

# Helper function
#' @export
set_epochs <- function(df, epoch = 30, subepochs = FALSE) {
    # Get quotient and remainder of euclidean division
    # of each time in seconds by 30.
    q <- df$Time %/% epoch
    df <- df %>% tibble::add_column(Epoch = as.factor(q + 1), .after = 1)
    if (subepochs) {
        r <- df$Time %% epoch
        qprime <- r %/% 5
        df <- df %>% tibble::add_column(Subepoch = as.factor(qprime + 1), .after = 2)
    }
    return(df)
}

#' Given an eeg object and a list of epochs, return a new EEG object 
#' with the said epochs removed.
#'
#' Notice that the extremes of the epoch are kept.
#'
#' @param object An eeg object.
#' @param epochs A vector of natural integers.
#' @return An EEG object.
#' @export
methods::setMethod(
    "drop_epochs",
    "eeg",
    function(object, epochs, epoch = 30) {
        df <- object@data %>% set_epochs(epoch)
        df <- droplevels(df[!df$Epoch %in% epochs, ])[-2]
        return(new("eeg", data = df, signals = object@signals))
    }
)

#' Given an eeg object, a list of epochs and
#' a list of subepochs, removes epoch-subepoch
#' pairs from the eeg data.
#' Note that elements in the epoch and subepoch
#' lists are assumed to have an element-wise
#' association. This means that, if the lists
#' are (e_1, ..., e_n) and (s_1, ..., s_n), then
#' epoch/subepoch pairs (e_i, s_i) are removed
#' for i in [1, n].
#'
#' @param object An eeg object.
#' @param epoch An natural number
#' @return An EEG object.
#' @export
methods::setMethod(
    "drop_subepochs",
    "eeg",
    function(object, epochs, subepochs, epoch = 30) {
        contaminated <- as.factor(paste(epochs, subepochs))
        df <- object@data %>% set_epochs(epoch, subepochs = TRUE)

        df <- df %>% tibble::add_column(
            Pairs = as.factor(paste(df$Epoch, df$Subepoch)),
            .after = "Subepoch"
        )
        df <- droplevels(df[!df$Pairs %in% contaminated, ])[-c(2, 3, 4)]
        return(new("eeg", data = df, signals = object@signals))
    }
)

#' Given an eeg object and an analysis object,
#' returns an artifact-rejected version of the
#' eeg.
#'
#' @param object An eeg object.
#' @param analysis An analysis object
#' @return An EEG object.
#' @export
methods::setMethod(
    "artf_reject",
    "eeg",
    function(object, analysis) {
        epoch_data <- extract_epochs(analysis)
        rejected <- drop_subepochs(object, epoch_data$Epoch, epoch_data$Subepoch)
        return(rejected)
    }
)
