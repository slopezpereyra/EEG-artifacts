library(tidyverse)

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
