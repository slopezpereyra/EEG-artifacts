#' @export
EEG <- R6::R6Class("EEG", list(
    data = tibble::tibble(),
    signals = tibble::tibble(),
    canoms = tibble::tibble(),
    panoms = tibble::tibble(),
    psd = tibble::tibble(),
    fs = 0,

    initialize = function(data_file, signals_file = NULL) {
        data <- readr::read_csv(data_file)
        if (!is.null(signals_file)) {
            signals <- readr::read_csv(signals_file)
            colnames(data)[-1] <- signals$Label %>% stringr::str_remove("EEG ") %>%
                stringr::str_remove("EOG")
        } else {
            signals <- tibble::tibble()
        }

        self$data <- data
        self$signals <- signals
        self$fs <- self$get_fs()
    },

    subset = function(s, e) {
        s_ind <- which(self$data$Time == s)
        e_ind <- which(self$data$Time == e)
        
        if (identical(s_ind, integer(0)) || identical(e_ind, integer(0))) {
            stop("You have not provided valid time bounds. Are you sure those values exist?")
        }
        self$data <- self$data[s_ind:e_ind, ]
    },

    resample = function(n) {
        self$data <- self$data[seq(1, nrow(self$data), n), ]
        self$fs = self$get_fs()
    },


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
    #' 
    get_epoch = function(which_epoch, epoch = 30) {
            self$subset(epoch * (which_epoch - 1), epoch * (which_epoch))
        },


    #' Returns the sampling frequency of the EEG object.
    #'
    #' @param object An eeg object.
    #'
    #' @return A number, the sampling frequency
    #' 
    get_fs = function() {
            delta_t <- self$data$Time[2] - self$data$Time[1]
            return(1 / delta_t)
        },
    
    #' Wrapper function for low-pass filtering
    #' a vector given a filtering frequency and
    #' a sampling frequency.
    #'
    #' 
    vlow_pass = function(vec, n, fs) {
        wpass <- n / (fs / 2) # Nyquist
        but <- gsignal::butter(5, wpass, "low", output = "Sos")
        low_pass <- gsignal::filter(but, unlist(vec))
        return(low_pass)
    },
    
    #' Wrapper function for high-pass filtering
    #' a vector given a filtering frequency and
    #' a sampling frequency.
    #'
    #' 
    vhigh_pass = function(vec, n, fs) {
        wpass <- n / (fs / 2) # Nyquist
        but <- gsignal::butter(5, wpass, "high", output = "Sos")
        high_pass <- gsignal::filter(but, unlist(vec))
        return(high_pass)
    },
    
    #' Wrapper function for bandpass filtering
    #' a vector given filtering frequencies and
    #' a sampling frequency.
    #'
    #' 
    vbandpass = function(vec, l, h, fs) {
        fpass <- c(l, h)
        wpass <- fpass / (fs / 2) # Nyquist
        but <- gsignal::butter(5, wpass, "pass", output = "Sos")
        pass <- gsignal::filter(but, unlist(vec))
        return(pass)
    },
    
    #' Given an EEG object and a numeric frequency n,
    #' applies a low-pass Butterworth filter to all channels.
    #'
    #' @param object An eeg object.
    #' @param n Filter frequency in Hz
    #'
    #' @return A new filtered EEG object
    #' 
    low_pass = function(n) {
            filt_df <- apply(self$data[-1],
                MARGIN = 2,
                FUN = function(x) self$vlow_pass(x, n, self$fs),
                simplify = FALSE
            ) %>%
                tibble::as_tibble() %>%
                tibble::add_column(Time = self$data$Time,
                                    .before = colnames(self$data)[2])
            self$data <- filt_df
        },
    
    
    #' Given an EEG object and a numeric frequency n,
    #' applies a high-pass Butterworth filter to all channels.
    #'
    #' @param object An eeg object.
    #' @param n Filter frequency in Hz
    #'
    #' @return A new filtered EEG object
    #' 
    high_pass = function(n) {
            filt_df <- apply(self$data[-1],
                MARGIN = 2,
                FUN = function(x) self$vhigh_pass(x, n, self$fs),
                simplify = FALSE
            ) %>%
                tibble::as_tibble() %>%
                tibble::add_column(Time = self$data$Time,
                                    .before = colnames(self$data)[2])
            self$data <- filt_df
        },
 

    #' Given an eeg object and a numeric frequency n,
    #' applies a bandpass filter to all channels.
    #'
    #' @param object An eeg object.
    #' @param l lower bound in Hz
    #' @param h higher bound in Hz
    #'
    #' @return A new filtered EEG object
    #' 
    bandpass = function(l, h) {
            filt_df <- apply(self$data[-1],
                MARGIN = 2,
                FUN = function(x) self$vbandpass(x, l, h, self$fs),
                simplify = FALSE
            ) %>%
                tibble::as_tibble() %>%
                tibble::add_column(Time = self$data$Time,
                                    .before = colnames(df)[2])
            self$data <- filt_df
        },
    
    
    #' Given an EEG object and a channel integer n,
    #' plots EEG record of the nth channel.
    #'
    #' @param object An eeg object.
    #' @param channel An integer indicating index of channel to plot.
    #'
    #' @return A ggplot object.
    #' 
    plot_channel = function(channel) {
            y <- self$data[-1][channel]
            p <- self$data %>%
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
                ggplot2::ylab(colnames(self$data[-1][channel]))
    
            return(p)
        },
    
    #' Given an EEG plot all the EEG channels
    #' in a vertical single-column layout.
    #'
    #' @param object An eeg object.
    #'
    #' @return A plot_grid object.
    #' 
    plot = function() {
            plots <- list()
            for (channel in 1:(ncol(self$data) - 1)) {
                p <- self$plot_channel(channel)
                plots[[channel]] <- p
            }
    
            return(cowplot::plot_grid(plotlist = plots, align = "v", ncol = 1))
        },
    
    #' Given an EEG, return its interactive visualization. This method should not
    #' be called on very large EEG objects (e.g., above an hour) without
    #' previous resampling. Interactive visualizations are computationally
    #' expensive.
    #'
    #' @param object An eeg object.
    #' @return A plotly figure.
    #' 
    iplot = function() {
            plots <- list()
            plots[[1]] <- plotly::plot_ly(
                self$data,
                type = "scatter",
                mode = "lines"
            ) %>%
                plotly::add_trace(
                    x = ~Time, y = unlist(self$data[, 2]),
                    name = colnames(self$data)[2]
                )
            for (i in 3:length(self$data)) {
                fig <- plotly::plot_ly(
                    self$data,
                    type = "scatter",
                    mode = "lines"
                ) %>%
                    plotly::add_trace(
                        x = ~Time, y = unlist(self$data[, i]),
                        name = colnames(self$data)[i]
                    )
                plots[[i - 1]] <- fig
            }
            return(plotly::subplot(
                plots,
                nrows = length(plots),
                shareX = TRUE
            ))
        },
    
    
    #' Given an eeg object and a list of epochs, return a new EEG object 
    #' with the said epochs removed.
    #'
    #' Notice that the extremes of the epoch are kept.
    #'
    #' @param object An eeg object.
    #' @param epochs A vector of natural integers.
    #' @return An EEG object.
    #' 
    drop_epochs = function(epochs, epoch = 30) {
            df <- self$data %>% set_epochs(epoch)
            df <- droplevels(df[!df$Epoch %in% epochs, ])[-2]
            self$data <- df
        },
    
    #' Given an eeg a list of epochs and
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
    #' 
    drop_subepochs = function(epochs, subepochs, epoch = 30) {
            contaminated <- as.factor(paste(epochs, subepochs))
            df <- self$data %>% set_epochs(epoch, subepochs = TRUE)
    
            df <- df %>% tibble::add_column(
                Pairs = as.factor(paste(df$Epoch, df$Subepoch)),
                .after = "Subepoch"
            )
            df <- droplevels(df[!df$Pairs %in% contaminated, ])[-c(2, 3, 4)]
            self$data <- df
        },

    # ---- Artifact related functions

    artf = function(alpha = 8, beta = 1) {
        print("Starting artifact analysis. This may take a couple of minutes...")
        start_time <- Sys.time()
        analysis <- anomaly::capa.mv(self$data[-1], type = "mean")
        canoms <- anomaly::collective_anomalies(analysis) %>%
            dplyr::filter(mean.change >= alpha) %>%
            set_timevars(data = self$data) %>%
            tibble::as_tibble()
        panoms <- anomaly::point_anomalies(analysis) %>%
            dplyr::filter(strength >= beta) %>%
            set_timevars(data = self$data) %>%
            tibble::as_tibble()
        end_time <- Sys.time()
        print(paste(
            "Analysis completed in ",
            (end_time - start_time), " seconds"
        ))
        self$canoms <- canoms
        self$panoms <- panoms
    },

    artf_stepwise = function(step_size = 30, alpha = 8) {
        print("Starting epoch-by-epoch artifact analysis. This may take a couple of minutes...")
        start_time <- Sys.time()
        # Set epochs for grouping
        t <- set_epochs(self$data, epoch = step_size) %>% head(-1)
        mps <- self$fs * step_size # measures per step
        grouped <- dplyr::group_by(t[-1], Epoch) %>%
            dplyr::group_map(~ anomaly::capa.mv(x = .x, type = "mean"))
        canoms <- grouped %>% lapply(function(x) anomaly::collective_anomalies(x) %>% dplyr::filter(mean.change >= alpha))
        panoms <- grouped %>% lapply(function(x) anomaly::point_anomalies(x))

        canoms <- mapply(function(x, y) x %>% dplyr::mutate(start = start + mps * (y - 1), end = end + mps * (y - 1)),
            canoms, seq_along(canoms),
            SIMPLIFY = FALSE
        ) %>%
            dplyr::bind_rows() %>%
            set_timevars(self$data)

        panoms <- mapply(function(x, y) x %>% mutate(location = location + mps * (y - 1)),
            panoms, seq_along(panoms),
            SIMPLIFY = FALSE
        ) %>%
            dplyr::bind_rows() %>%
            set_timevars(self$data)
        
        end_time <- Sys.time()
        print(paste(
            "Analysis completed in ",
            (end_time - start_time), " seconds"
        ))
        self$canoms <- canoms
        self$panoms <- panoms
    },

    get_contaminated_channels = function() {
        chans <- union(
            unique(self$canoms$variate),
            unique(self$panoms$variate)
        )
        return(chans)
    },

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
    #' 
    set_plot_data = function(chan, data) {
            canoms <- dplyr::filter(self$canoms, variate == chan)
            panoms <- dplyr::filter(self$panoms, variate == chan)
            # Get all indexes between start and end of canoms
            locations <- mapply(function(x, y) x:y, canoms$start, canoms$end)
            # Unite with point anomalies
            locations <- union(unlist(locations), unlist(panoms$location)) %>%
                        as.integer()
            time_of_anomalies <- lubridate::as_datetime(unlist(data[locations, 1]))
            values <- unlist(data[locations, chan + 1])
            df <- tibble::tibble(A = time_of_anomalies, B = values)
            return(df)
        },
    
    #' Plots analysis results for a specific channel.
    #'
    #' @param eeg An Analysis object.
    #' @param chan An integer that points to the EEG channel.
    #' @param size (optional) The size of the red dots that signal an anomaly.
    #' @return A data frame
    #' 
    plot_channel_artifacts = function(chan, size = 0.2) {
            df <- self$set_plot_data(chan, self$data)
            eeg <- self$plot_channel(channel = chan)
            p <- eeg +
                ggplot2::geom_point(
                    data = df, ggplot2::aes(A, B),
                    inherit.aes = FALSE, color = "red",
                    size = size
                )
            return(p)
        },

    #' Plots all channels with any kind of anomalies.
    #'
    #' @param x An Analysis eeg.
    #' @return A plot_grid eeg.
    #' 
    plot_artifacts = function(size = 0.2) {
            plots <- list()
            channels <- self$get_contaminated_channels()
            for (channel in channels) {
                p <- self$plot_channel_artifacts(channel, size)
                plots[[channel]] <- p
            }
            return(cowplot::plot_grid(plotlist = plots, align = "v", ncol = 1))
        },


    #' Filters an analysis eeg so as to keep only those
    #' anomalies whose normalized strength is greater than x.
    #'
    #' @param eeg An analysis object
    #' @param x The strength threshold
    #' @param f A function used to normalize strengths (otherwise collective
    #' and point anomalies are incomparable). Defaults to min-max normalization.
    #'
    #' @return An analysis eeg.
    #' 
    sfilter = function(x, f = minmax_normalization) {
            self$canoms <- f(self$canoms) %>% dplyr::filter(mean.change >= x)
            self$panoms <- f(self$panoms) %>% dplyr::filter(strength >= x)
        },

    #' Helper function
    #' 

    #' Given an analysis eeg, returns a data frame with the epoch-subepoch pairs
    #' found to contain anomalies, and their respective average strength.
    #'
    #' @param eeg An analysis object.
    #' @return A data frame containing epoch-subepoch pairs
    #' found to contain artifacts as well as the average
    #' strength of the artifacts for each subepoch.
    #' 
    extract_epochs = function() {
            l <- list(self$canoms, self$panoms)
            i <- which(lapply(l, nrow) == 0)
            # If both data frames have length 0
            if (length(i) > 1) {
                stop("Empty analysis error")
            }
            # If both data frames are non-empty
            if (identical(i, integer(0))) {
                a <- canoms_avg_epoch_strength(self$canoms)
                b <- panoms_avg_epoch_strength(self$panoms)
                return(dplyr::full_join(a, b))
            }
            # If only one is non-empty
            if (i == 1) {
                return(panoms_avg_epoch_strength(self$panoms))
            }
            return(canoms_avg_epoch_strength(self$canoms))
        },

    #' Given an eeg object and an analysis object,
    #' returns an artifact-rejected version of the
    #' eeg.
    #'
    #' @param object An eeg object.
    #' @param analysis An analysis object
    #' @return An EEG object.
    #' 
    artf_reject = function() {
            epoch_data <- self$extract_epochs()
            clone <- self$clone()
            clone$drop_subepochs(epoch_data$Epoch, epoch_data$Subepoch)
            clone$canoms <- clone$panoms <- tibble::tibble()
            return(clone)
        },

    # --- PSD ---


    spectrogram = function(channel, max_freq = 30, freq = 4) {
        rsleep::spectrogram(unlist(self$data[-1][channel]),
                            sRate = self$fs, maxFreq = max_freq,
                            freq = freq)
    },

    channel_psd = function(channel) {
        welch <- gsignal::pwelch(as.matrix(self$data[-1][channel]), fs = self$fs)
        psd <- welch$spec %>%
            apply(log10, MARGIN = 2) %>%
            tibble::as_tibble()
        psd$Fqc <- welch$freq
        return(psd)
    },

    compute_psd = function() {
        pwelch <- gsignal::pwelch(as.matrix(self$data[-1]), fs = self$fs)
        psd <- pwelch$spec %>%
            apply(log10, MARGIN = 2) %>%
            tibble::as_tibble()
        psd$Fqc <- pwelch$freq
        self$psd <- psd
    },
    #' Given an PSD data frame as returned by the psd(eeg, ...) function,
    #' plot the spectrum of all channels.
    #'
    #' @param object An data frame as returned by the psd(eeg, ...) function.
    #' @return A ggplot object.
    #'
    #' 
    plot_psd = function(xlim = 250) {
        tall <- reshape2::melt(self$psd, id.vars = "Fqc")
        p <- ggplot2::ggplot(tall, ggplot2::aes(Fqc, value, col=variable)) +
            ggplot2::geom_line() +
            ggplot2::xlim(c(0, xlim))
        return(p)
    },

    iplot_psd = function(xlim = 250) {
        psd <- reshape2::melt(self$psd, id.vars = "Fqc")
        fig <- plotly::plot_ly(
            psd,
            type = "scatter",
            mode = "lines"
        ) %>%
            plotly::add_trace(
                x = ~Fqc, y = ~value, color = ~variable
            ) %>%
            layout(
                xaxis = list(
                    title = "Frequency in Hz",
                    zeroline = F,
                    range = c(0, xlim)
                ),
                yaxis = list(
                    title = "log10 Spectrum",
                    zeroline = F
                )
            )
        return(fig)
    }
    )
)
