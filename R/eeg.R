library(gsignal)
library(rsleep)

#' EEG class.
#' @slot data A data frame containing EEG records
#' @slot signals A data frame containing signal information (may be empty).
#' @export
setClass("eeg",
    slots = list(
        data = "data.frame",
        signals = "data.frame"
    )
)

#' Show method for the EEG class that prints the data slot.
#' @param object An EEG object.
#' @export
setMethod(
    "show",
    "eeg",
    function(object) {
        print(object@data)
    }
)

# ---------- EEG Generics ------------

#' @export
setGeneric("plot_channel", function(object, channel) standardGeneric("plot_channel"))

#' @export
setGeneric("bandpass", function(object, l, h) standardGeneric("bandpass"))

#' @export
setGeneric("high_pass", function(object, n) standardGeneric("high_pass"))

#' @export
setGeneric("subset_eeg", function(object, s, e) standardGeneric("subset_eeg"))

#' @export
setGeneric("resample_eeg", function(object, n) standardGeneric("resample_eeg"))

#' @export
setGeneric("get_epoch", function(object, which_epoch, epoch = 30) standardGeneric("get_epoch"))

#' @export
setGeneric("get_samples_in_epoch", function(object, epoch) standardGeneric("get_samples_in_epoch"))

#' @export
setGeneric("get_sampling_frequency", function(object) standardGeneric("get_sampling_frequency"))

#' @export
setGeneric("low_pass", function(object, n) standardGeneric("low_pass"))

#' @export
setGeneric("remove_epoch", function(object, epoch) standardGeneric("remove_epoch"))

#' Read a .csv data file containing EEG data and an optionall
#' signals file and return an eeg object. If the signals file is provided,
#'  channel names are appropriately set.
#'
#' @param data_file .csv file containing eeg data.
#' @param signals_file .csv file containing signal information
#'
#' @return An eeg object.
#' @export
load_eeg <- function(data_file, signals_file = NULL) {
    data <- read_csv(data_file)
    if (!is.null(signals_file)) {
        signals <- read_csv(signals_file)
        colnames(data)[-1] <- signals$Label %>%
            str_remove("EEG ") %>%
            str_remove("EOG")
    } else {
        signals <- tibble()
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
setMethod(
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



#' Given an eeg object and an integer n,
#' subset the EEG's data by keeping only one every
#' n values and return a new EEG with the subsetted data.
#'
#' Logically, the new EEG will have a length of 1/n
#' the length of the subsetted EEG.
#'
#' @param object An eeg object.
#' @param n An integer.
#'
#' @return An EEG object
#' @export
setMethod(
    "resample_eeg",
    "eeg",
    function(object, n) {
        df <- object@data[seq(1, nrow(object@data), n), ]
        return(new("eeg", data = df, signals = object@signals))
    }
)
#' Returns new eeg object whose data is a subset of given eeg
#' on the given epoch.
#'
#' This function wraps a specific call of subset_eeg.
#'
#' @param object An eeg object.
#' @param which_epoch What epoch to subset.
#' @param epoch How many seconds make up an epoch? Defaults to 30.
#'
#' @return A new eeg whose data is the requested epoch.
#' @export
setMethod(
    "get_epoch",
    "eeg",
    function(object, which_epoch, epoch = 30) {
        subset_eeg(object, epoch * (which_epoch - 1), epoch * (which_epoch))
    }
)


#' Given an eeg object, determine the number of values
#' that make up epoch seconds.
#'
#' @param object An eeg object.
#' @param epoch A time in seconds.
#'
#' @return An integer representing the number of values
#' that make up a time-frame of length epoch.
#' @export
setMethod(
    "get_samples_in_epoch",
    "eeg",
    function(object, epoch) {
        return(which(object@data["Time"] == epoch))
    }
)

#' Returns sampling frequency of the EEG object
#'
#' @param object An eeg object.
#'
#' @return An number, the sampling frequency
#' @export
setMethod(
    "get_sampling_frequency",
    "eeg",
    function(object) {
        delta_t <- object@data$Time[2] - object@data$Time[1]
        return(1 / delta_t)
    }
)


#' Given an eeg object and a numeric frequency n,
#' applies a low-pass Butterworth filter.
#'
#' @param object An eeg object.
#' @param n Filter frequency in Hz
#'
#' @return A new filtered EEG object
#' @export
setMethod(
    "low_pass",
    "eeg",
    function(object, n) {
        for (chan in 1:(ncol(object@data) - 1)) {
            fs <- get_sampling_frequency(object)
            fs <- get_sampling_frequency(object)
            fpass <- n
            wpass <- fpass / (fs / 2) # Nyquist
            but <- gsignal::butter(3, wpass, "low")
            highp <- gsignal::filter(but, unlist(object@data[-1][chan]))
            object@data[-1][chan] <- highp
        }
        return(object)
    }
)



#' Given an eeg object and a numeric frequency n,
#' applies high pass Butterworth filter.
#'
#' @param object An eeg object.
#' @param n Filter frequency in Hz.
#'
#' @return A new filtered EEG object
#' @export
setMethod(
    "high_pass",
    "eeg",
    function(object, n) {
        for (chan in 1:(ncol(object@data) - 1)) {
            fs <- get_sampling_frequency(object)
            fs <- get_sampling_frequency(object)
            fpass <- n
            wpass <- fpass / (fs / 2) # Nyquist
            but <- gsignal::butter(3, wpass, "high")
            highp <- gsignal::filter(but, unlist(object@data[-1][chan]))
            object@data[-1][chan] <- highp
        }
        return(object)
    }
)


#' Given an eeg object and a numeric frequency n,
#' applies a bandpass filter.
#'
#' @param object An eeg object.
#' @param l lower bound in Hz
#' @param h higher bound in Hz
#'
#' @return A new filtered EEG object
#' @export
setMethod(
    "bandpass",
    "eeg",
    function(object, l, h) {
        for (chan in 1:(ncol(object@data) - 1)) {
            fs <- get_sampling_frequency(object)
            fpass <- c(l, h)
            wpass <- fpass / (fs / 2) # Nyquist
            but <- gsignal::butter(5, wpass, "pass")
            pass <- gsignal::filter(but, unlist(object@data[-1][chan]))
            object@data[-1][chan] <- pass
        }
        return(object)
    }
)


#' Given an eeg object and a channel integer n,
#' plots EEG record of the nth channel.
#'
#' @param object An eeg object.
#' @param channel An integer indicating index of channel to plot.
#'
#' @return A ggplot object.
#' @export
setMethod(
    "plot_channel",
    "eeg",
    function(object, channel) {
        y <- object@data[-1][channel]
        p <- object@data %>%
            mutate(Time = as_datetime(Time)) %>%
            ggplot(
                aes(
                    Time,
                    unlist(y)
                )
            ) +
            geom_line() +
            scale_x_datetime(date_labels = "%H:%M:%S") +
            xlab("") +
            ylab(colnames(object@data[-1][channel]))

        return(p)
    }
)

#' Given an eeg object, plot all the EEG channels
#' in a vertical single-column layout.
#'
#' @param object An eeg object.
#'
#' @return A plot_grid object.
#' @export
setMethod(
    "plot",
    "eeg",
    function(x) {
        plots <- list()
        for (channel in 1:(ncol(x@data) - 1)) {
            p <- plot_channel(x, channel)
            plots[[channel]] <- p
        }

        return(plot_grid(plotlist = plots, align = "v", ncol = 1))
    }
)

#' Given an eeg object and an epoch, removes
#' all samples from that epoch from the EEG data.
#'
#' Notice that the extremes of the epoch are kept.
#'
#' @param object An eeg object.
#' @param epoch An natural number
#' @return An EEG object.
#' @export
setMethod(
    "remove_epoch",
    "eeg",
    function(object, epoch) {
        s_ind <- which(object@data$Time == 30 * (epoch - 1)) + 1
        e_ind <- which(object@data$Time == 30 * epoch) - 1
        df <- object@data[-c(s_ind:e_ind), ]
        return(new("eeg", data = df, signals = object@signals))
    }
)
