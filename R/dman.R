
#' Create an empty data frame to be filled with epoch-subepoch
#' anomalous pairs during stepwise analysis.
#'
#' @return An empty dataframe.
#' @export
create_epoch_data <- function() {
    results <- tibble::tibble(
        Epoch = numeric(),
        Subepoch = numeric()
        # channels = list(),
        # segment_strength = numeric(),
        # point_strength = numeric()
    )
    return(results)
}


#' Update a data frame containing anomalous epoch-subepoch
#' pairs given a certain analysis results.
#'
#' @return A data frame as defined by create_epoch_data().
#' @return An analysis object.
#' @export
update_epochs <- function(epoch_data, analysis) {
    canoms <- analysis@canoms
    panoms <- analysis@panoms

    collective_epochs <- unique(paste(canoms$Epoch, canoms$Subepoch))
    point_epochs <- unique(paste(panoms$Epoch, panoms$Subepoch))

    epochs <- union(collective_epochs, point_epochs)
    epoch_list <- lapply(strsplit(epochs, " "), as.numeric)
    # chans <- union(unique(canoms$variate), unique(panoms$variate))

    epoch_data <- tibble::add_row(epoch_data,
        Epoch = unlist(lapply(epoch_list, `[[`, 1)),
        Subepoch = unlist(lapply(epoch_list, `[[`, 2))
    )

    return(epoch_data)
}


#' @export
methods::setGeneric("extract_epochs", function(object) standardGeneric("extract_epochs"))

#' data attribute is such subset.
#'
#' @param object An eeg object.
#' @param s Starting time of the subset in seconds.
#' @param e Ending time of the subset in seconds.
#'
#' @return A new eeg whose data is the subset ranging from
#' second s to e of the object's data attribute.
#' @export
methods::setMethod(
    "extract_epochs",
    "analysis",
    function(object) {
        df <- create_epoch_data()
        df <- update_epochs(df, object)
        return(df)
    }
)
