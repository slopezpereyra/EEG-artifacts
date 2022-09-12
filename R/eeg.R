library(signal)

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

#' @export
setGeneric(
    "partition.eeg",
    function(object, s, e) {
        standardGeneric("partition.eeg")
    }
)

#' Given an eeg object, determine the subset of
#' its data attribute that contains records from
#' second s to e, and return a new eeg object whose
#' data attribute is such subset.
#'
#' @param object An eeg object.
#' @param s Starting time of the subset in seconds.
#' @param e Ending time of the subset in seconds.
#'
#' @return A new eeg whose data is the subset ranging from
#' second s to e of the object's data attribute.
#' @export
setMethod(
    "partition.eeg",
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

#' @export
setGeneric(
    "low.pass",
    function(object, n) {
        standardGeneric("low.pass")
    }
)

#' Given an eeg object and a numeric frequency n,
#' applies a 1/n Hz low-pass Butterworth filter.
#'
#' @param object An eeg object.
#' @param n Filter frequency.
#' @import signal
#'
#' @return A new filtered EEG object
#' @export
setMethod(
    "low.pass",
    "eeg",
    function(object, n) {
        for (chan in 1:(ncol(object@data) - 1)) {
            bf <- butter(3, 1 / n) # 1/n Hz low-pass filter
            lowp <- signal::filter(bf, unlist(object@data[-1][chan]))
            object@data[-1][chan] <- lowp
        }
        return(object)
    }
)

#' @export
setGeneric(
    "lower.res",
    function(object, n) {
        standardGeneric("lower.res")
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
#' @return A new filtered EEG object
#' @export
setMethod(
    "lower.res",
    "eeg",
    function(object, n) {
        i <- c(1) # Include one since following computation skips it.
        if (n == 1) {
            return(object)
        }
        for (x in seq_len(nrow(object@data) / n)) {
            i <- append(i, x * n + 1)
        }
        lowered <- object@data[i, ]
        return(new("eeg", data = lowered, signals = object@signals))
    }
)

#' @export
setGeneric(
    "plot.channel",
    function(object, channel) {
        standardGeneric("plot.channel")
    }
)

#' Given an eeg object and a channel integer n,
#' plots EEG record of the nth channel.
#'
#' @param object An eeg object.
#' @param channel An integer indicating index of channel to plot.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @export
setMethod(
    "plot.channel",
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

#' @export
setGeneric(
    "plot.eeg",
    function(object) {
        standardGeneric("plot.eeg")
    }
)

#' Given an eeg object, plot all the EEG channels
#' in a vertical single-column layout.
#'
#' @param object An eeg object.
#'
#' @return A plot_grid object.
#' @importFrom cowplot plot_grid
#' @export
setMethod(
    "plot.eeg",
    "eeg",
    function(object) {
        plots <- list()
        for (channel in 1:(ncol(object@data) - 1)) {
            p <- plot.channel(object, channel)
            plots[[channel]] <- p
        }

        return(plot_grid(plotlist = plots, align = "v", ncol = 1))
    }
)

#' Read a .csv data file containing EEG data and an optionall
#' signals file and return an eeg object. If the signals file is provided,
#'  channel names are appropriately set.
#'
#' @param data_file .csv file containing eeg data.
#' @param signals_file .csv file containing signal information
#'
#' @return An eeg object.
#' @export
load.eeg <- function(data_file, signals_file = NULL) {
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


#' Create an empty data frame to be filled with epoch-subepoch
#' anomalous pairs during stepwise analysis.
#'
#' @return An empty dataframe.
create.epoch.data <- function() {
    results <- tibble(
        epoch = numeric(),
        subepoch = numeric()
        # channels = list(),
        # segment_strength = numeric(),
        # point_strength = numeric()
    )
    return(results)
}


#' Update a data frame containing anomalous epoch-subepoch
#' pairs given a certain analysis results.
#'
#' @return A data frame as defined by create.epoch.data().
#' @return An analysis object.
update.epochs <- function(epoch_data, analysis) {
    canoms <- analysis@canoms
    panoms <- analysis@panoms

    collective_epochs <- unique(paste(canoms$Epoch, canoms$Subepoch))
    point_epochs <- unique(paste(panoms$Epoch, panoms$Subepoch))

    epochs <- union(collective_epochs, point_epochs)
    epoch_list <- lapply(strsplit(epochs, " "), as.numeric)
    chans <- union(unique(canoms$variate), unique(panoms$variate))

    epoch_data <- add_row(epoch_data,
        epoch = unlist(lapply(epoch_list, `[[`, 1)),
        subepoch = unlist(lapply(epoch_list, `[[`, 2))
    )

    return(epoch_data)
}
