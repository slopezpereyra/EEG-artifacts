

# Scripts on this file are dedicated to plotting EEG data
# and M-CAPA analysis results.

library(ggplot2)
library(cowplot)
library(lubridate)


setClass("analysis",
    slots = list(
        canoms = "data.frame",
        panoms = "data.frame",
        origin = "data.frame"
    )
)

setMethod(
    "show",
    "analysis",
    function(object) {
        View(object@canoms)
        View(object@panoms)
        View(object@origin)
    }
)

setGeneric(
    "has.anomalies",
    function(object, step, ...) {
        standardGeneric("has.anomalies")
    }
)

# " Check if analysis as returned by the analyze() function
# " found anomalies or not.
# " @param analysis An analysis as returned by the analyze() function.

setMethod(
    "has.anomalies",
    "analysis",
    function(object) {
        return(nrow(object@canoms) > 0 | nrow(object@panoms) > 0)
    }
)


setGeneric(
    "get.anomalous.channels",
    function(object, channel) {
        standardGeneric("get.anomalous.channels")
    }
)

setMethod(
    "get.anomalous.channels",
    "analysis",
    function(object) {
        chans <- union(
            unique(object@canoms$variate),
            unique(object@panoms$variate)
        )
        return(chans)
    }
)


setGeneric(
    "draw.eeg",
    function(object, channel) {
        standardGeneric("draw.eeg")
    }
)

setMethod(
    "draw.eeg",
    "analysis",
    function(object, channel) {
        plot <- object@origin %>%
            mutate(Time = as_datetime(Time)) %>%
            ggplot(
                aes(
                    Time,
                    unlist(object@origin[-1][channel])
                )
            ) +
            geom_line() +
            scale_x_datetime(date_labels = "%H:%M:%S") +
            xlab("") +
            ylab(colnames(object@origin[-1][channel]))
    }
)

setGeneric(
    "draw.clusters",
    function(object, channel) {
        standardGeneric("draw.clusters")
    }
)

setMethod(
    "draw.clusters",
    "analysis",
    function(object, channel) {
        canoms <- object@canoms %>%
            dplyr::filter(variate == channel)
        # mutate(
        #    start = as_datetime(object@origin$Time[object@canoms$start])
        # )
        areas <- tibble(
            xmin = as_datetime(canoms$Time),
            xmax = as_datetime(object@origin$Time[canoms$end]),
            ymin = -315, ymax = 315, alpha = 0.3
        )
        print(canoms)
        print(areas)
        eeg <- draw.eeg(object, channel = channel)
        plot <- eeg +
            geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = 0.3, fill = "red",
                data = transform(areas, as.character(seq_len(nrow(areas)))),
                inherit.aes = FALSE
            ) + xlab("") + ylab(colnames(object@origin[-1][channel]))
        return(plot)
    }
)

setGeneric(
    "draw.points",
    function(object, channel) {
        standardGeneric("draw.points")
    }
)

setMethod(
    "draw.points",
    "analysis",
    function(object, channel) {
        points <- object@panoms %>%
            dplyr::filter(variate == channel)
        point_values <- unlist(object@origin[-1][points$location, channel])
        points <- tibble(A = as_datetime(points$Time), B = point_values)

        eeg <- draw.eeg(object, channel = channel)
        plot <- eeg +
            geom_point(data = points, aes(A, B), inherit.aes = FALSE, color = "red") +
            xlab("") + ylab(colnames(object@origin[-1][channel]))

        return(plot)
    }
)

setGeneric(
    "draw.anomalies",
    function(object, channel) {
        standardGeneric("draw.anomalies")
    }
)

setMethod(
    "draw.anomalies",
    "analysis",
    function(object, channel) {
        points <- object@panoms %>%
            dplyr::filter(variate == channel)
        point_values <- unlist(object@origin[-1][points$location, channel])
        points <- tibble(A = as_datetime(points$Time), B = point_values)

        canoms <- object@canoms %>% dplyr::filter(variate == channel)
        areas <- tibble(
            xmin = as_datetime(canoms$Time),
            xmax = as_datetime(object@origin$Time[canoms$end]),
            ymin = -300, ymax = 300, alpha = 0.3
        )

        eeg <- draw.eeg(object, channel = channel)
        plot <- eeg +
            geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                alpha = 0.3, fill = "red",
                data = transform(areas, as.character(1:nrow(areas))),
                inherit.aes = FALSE
            ) +
            geom_point(data = points, aes(A, B), inherit.aes = FALSE, color = "red") +
            xlab("") + ylab(colnames(object@origin[-1][channel]))

        return(plot)
    }
)

setGeneric(
    "plot.channel",
    function(object, channel) {
        standardGeneric("plot.channel")
    }
)

setMethod(
    "plot.channel",
    "analysis",
    function(object, channel) {
        canoms <- object@canoms %>% dplyr::filter(variate == channel)
        point_locs <- object@panoms %>%
            dplyr::filter(variate == channel) %>%
            .$location
        if (nrow(canoms) == 0 & length(point_locs) == 0) {
            return(draw.egg(object, channel))
        } else if (nrow(canoms) == 0) {
            return(draw.points(object, channel))
        } else {
            return(draw.anomalies(object, channel))
        }
        return(plot)
    }
)

setGeneric(
    "plot.channels",
    function(object) {
        standardGeneric("plot.channels")
    }
)

setMethod(
    "plot.channels",
    "analysis",
    function(object) {
        plots <- list()
        channels <- get.anomalous.channels(object)
        for (channel in channels)
        {
            p <- plot.channel(object, channel)
            plots[[channel]] <- p
        }
        return(plot_grid(plotlist = plots, align = "v", ncol = 1))
    }
)

setGeneric(
    "plot",
    function(object, save = FALSE) {
        standardGeneric("plot")
    }
)

setMethod(
    "plot",
    "analysis",
    function(object, save = FALSE) {
        anom_plot <- plot.channels(object)
        if (save == TRUE) {
            start <- seconds_to_period(object@origin$Time[1]) %>% format.time()
            end <- seconds_to_period(tail(object@origin$Time, 1)) %>% format.time()
            time <- paste(start, end, sep = " to ")
            ggsave(paste("/home/santi/work/EEG-artifacts/images/", time, ".png", sep = ""), plot = anom_plot)
        }
        return(anom_plot)
    }
)
