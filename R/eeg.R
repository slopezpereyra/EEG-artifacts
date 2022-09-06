library(tidyverse)
library(signal)

setClass("eeg",
    slots = list(
        data = "data.frame",
        signals = "data.frame"
    )
)

setMethod(
    "show",
    "eeg",
    function(object) {
        print(object@data)
    }
)

setMethod(
    "na.omit",
    "eeg",
    function(object) {
        df <- na.omit(object@data)
        return(new("eeg", data = df, signals = object@signals))
    }
)

setGeneric(
    "partition.eeg",
    function(object, s, e) {
        standardGeneric("partition.eeg")
    }
)

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

setGeneric(
    "low.pass",
    function(object, n) {
        standardGeneric("low.pass")
    }
)

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

setGeneric(
    "lower.res",
    function(object, n) {
        standardGeneric("lower.res")
    }
)

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

setGeneric(
    "draw.channel",
    function(object, channel) {
        standardGeneric("draw.channel")
    }
)

setMethod(
    "draw.channel",
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

setGeneric(
    "draw",
    function(object) {
        standardGeneric("draw")
    }
)

setMethod(
    "draw",
    "eeg",
    function(object) {
        plots <- list()
        for (channel in 1:(ncol(object@data) - 1)) {
            p <- draw.channel(object, channel)
            plots[[channel]] <- p
        }
        return(plot_grid(plotlist = plots, align = "v", ncol = 1))
    }
)

load.eeg <- function(data_file, signals_file) {
    data <- read_csv(data_file)
    signals <- read_csv(signals_file)
    colnames(data)[-1] <- signals$Label %>%
        str_remove("EEG ") %>%
        str_remove("EOG")

    return(new("eeg", data = data, signals = signals))
}

create.epoch.data <- function() {
    results <- tibble(
        epoch = numeric(),
        subepoch = numeric(),
        channels = list(),
        segment_strength = numeric(),
        point_strength = numeric()
    )
    return(results)
}

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
