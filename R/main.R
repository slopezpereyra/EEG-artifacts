# EEG R6 class

#' @title EEG
#'
#' @description
#' Abstract base class for EEG data.
#' @export
EEG <- R6::R6Class("EEG", public = list(
    #' @field data (`tibble`)\cr
    #' Data frame (tibble) with the EEG data.
    data = tibble::tibble(),
    #' @field canoms (`tibble`)\cr
    #' Data frame (tibble) with the collective anomalies data.
    canoms = tibble::tibble(),
    #' @field panoms (`tibble`)\cr
    #' Data frame (tibble) with the point anomalies data.
    panoms = tibble::tibble(),
    #' @field psd (`tibble`)\cr
    #' Data frame (tibble) with the power spectral density data.
    psd = tibble::tibble(),
    #' @field spindles (`tibble`)\cr
    #' Data frame (tibble) with the spindle data.
    spindles = tibble::tibble(),
    #' @field fs (`numeric`)
    #'  Sampling frequency of the EEG.
    fs = 0,

    #' @description
    #' Creates a new instance of the EEG class.
    #'
    #' @param data_file (`string`)
    #'   Directory of .csv file containing EEG data.
    #' @param set_epochs (`bool = TRUE`)
    #'   Whether to set the epoch and subepoch features upon reading the 
    #'   EEG data. Defaults to TRUE, but should be FALSE if the EEG data 
    #'   being read already has these features.
    #' @param epoch (`int`)
    #'   How many seconds make up an epoch? Only relevant if set_epochs is TRUE.
    initialize = function(data_file, 
                          set_epochs=TRUE, 
                          epoch = 30) {
        f <- ifelse(endsWith(data_file, ".edf"), read_edf, readr::read_csv)
        self$data <- f(data_file)
        if (set_epochs) {
            self$data <- set_epochs(self$data, epoch, subepochs = TRUE)
        }
        self$fs <- self$get_fs()
        gc() # ::edf module seems to produce memory ovefload
    },

    #' @description
    #' Subsets EEG data to keep only measures in the interval [s, e], where
    #' `s` and `e` are times in seconds and `e > s`.
    #'
    #' @param s (`integer`).
    #' @param e (`integer`).
    #' @return void
    subset_by_seconds = function(s, e) {
        s_ind <- which(self$data$Time == s)
        e_ind <- which(self$data$Time == e)

        if (identical(s_ind, integer(0)) || identical(e_ind, integer(0))) {
            stop("You have not provided valid time bounds. Are you sure those
                 values exist?")
        }
        self$data <- self$data[s_ind:e_ind, ]
    },

    #' @description
    #' Subsets EEG data to keep only measures in the interval [s, e], where
    #' `s` and `e` are epochs and `e >= s`. If `e == s` then only data from
    #' that single epoch is kept.
    #'
    #' @param s (`integer`).
    #' @param e (`integer`).
    #' @return void
    subset = function(s, e) {
        self$data <- self$data %>%
            dplyr::filter(Epoch %in% c(s:e))
    },

    #' @description
    #' Subsets EEG data to keep only one for every `n` samples. This is a
    #' brute resampling algorithm and should only not be used except for
    #' special purposes.
    #'
    #' @param n (`integer`).
    #' @return void
    resample = function(n) {
        self$data <- self$data[seq(1, nrow(self$data), n), ]
        self$fs <- self$get_fs()
    },

    #' @description
    #' Returns the sampling frequency of the EEG data. The frequency is computed
    #' from the data and not extracted from the $fs field. This method is useful
    #' for updating the $fs field after resampling or other data modifications.
    #'
    #' @return void
    get_fs = function() {
        delta_t <- self$data$Time[2] - self$data$Time[1]
        return(1 / delta_t)
    },

    #' @description
    #' Private (do not use)
    #'
    #' @param vec (`vector`).
    #' @param n (`int`).
    #' @param fs (`int`).
    #'
    #' @return void
    vlow_pass = function(vec, n, fs) {
        wpass <- n / (fs / 2) # Nyquist
        but <- gsignal::butter(5, wpass, "low", output = "Sos")
        low_pass <- gsignal::filter(but, unlist(vec))
        return(low_pass)
    },


    #' @description
    #' Private (do not use)
    #'
    #' @param vec (`vector`).
    #' @param n (`int`).
    #' @param fs (`int`).
    #'
    #' @return void
    vhigh_pass = function(vec, n, fs) {
        wpass <- n / (fs / 2) # Nyquist
        but <- gsignal::butter(5, wpass, "high", output = "Sos")
        high_pass <- gsignal::filter(but, unlist(vec))
        return(high_pass)
    },

    #' @description
    #' Private (do not use)
    #'
    #' @param vec (`vector`).
    #' @param n (`int`).
    #' @param l (`int`).
    #' @param h (`int`).
    #' @param fs (`int`).
    #'
    #' @return void
    vbandpass = function(vec, l, h, fs) {
        fpass <- c(l, h)
        wpass <- fpass / (fs / 2) # Nyquist
        but <- gsignal::butter(5, wpass, "pass", output = "Sos")
        pass <- gsignal::filter(but, unlist(vec))
        return(pass)
    },

    #' @description
    #' Applies a low-pass filter to the EEG data. The value of `n` is used
    #' to compute the pass filter `F = n / (fs / 2)` (Nyquist frequency).
    #'
    #' @param n (`numeric`).
    #' @return void
    low_pass = function(n) {
        filt_df <- apply(self$data[, -c(1:3)],
            MARGIN = 2,
            FUN = function(x) self$vlow_pass(x, n, self$fs),
            simplify = FALSE
        ) %>%
            tibble::as_tibble() %>%
            tibble::add_column(Time = self$data$Time,
                                .before = colnames(self$data)[2])
        self$data <- filt_df
    },


    #' @description
    #' Applies a high-pass filter to the EEG data. The value of `n` is used
    #' to compute the pass filter `F = n / (fs / 2)` (Nyquist frequency).
    #'
    #' @param n (`numeric`).
    #' @return void
    high_pass = function(n) {
        filt_df <- apply(self$data[, -c(1:3)],
            MARGIN = 2,
            FUN = function(x) self$vhigh_pass(x, n, self$fs),
            simplify = FALSE
        ) %>%
            tibble::as_tibble() %>%
            tibble::add_column(Time = self$data$Time,
                                .before = colnames(self$data)[2])
        self$data <- filt_df
    },


    #' @description
    #' Applies a bandpass filter to the EEG data. Parameters `l` and `h` are
    # ' the low and high frequency bounds.
    #'
    #' @param l (`numeric`).
    #' @param h (`numeric`).
    #' @return void
    bandpass = function(l, h) {
        filt_df <- apply(self$data[, -c(1:3)],
            MARGIN = 2,
            FUN = function(x) self$vbandpass(x, l, h, self$fs),
            simplify = FALSE
        ) %>%
            tibble::as_tibble() %>%
            tibble::add_column(Time = self$data$Time,
                                .before = colnames(df)[2])
        self$data <- filt_df
    },


    #' @description
    #' Plots the EEG channel number `channel`. If `s` and `e` are non-zero the
    #' plot only shows EEG data from epoch `s` to epoch `e` (inclusive). If
    #  both are zero (default) the full EEG channel is ploted.
    #'
    #' @param channel (`integer`).
    #' @param s (`integer`) = 0
    #' @param e (`integer`) = 0
    #' @return ggplot
    plot_channel = function(channel, s = 0, e = 0) {
        data <- self$data
        if (s != 0 && e != 0) {
            data <- dplyr::filter(data, Epoch %in% c(s:e))
        }
        y <- data[, -c(1:3)][channel]
        p <- data %>%
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
            ggplot2::ylab(colnames(data[, -c(1:3)][channel]))

        return(p)
    },

    #' @description
    #' Plots the EEG data.
    #'
    #' @return cowplot
    plot = function() {
        plots <- list()
        for (channel in 1:(ncol(self$data) - 3)) {
            p <- self$plot_channel(channel)
            plots[[channel]] <- p
        }

        return(cowplot::plot_grid(plotlist = plots, align = "v", ncol = 1))
    },

    #' @description
    #' Produces an interactive PlotLy plot of the EEG data.
    #'
    #' @return plotly
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


    #' @description
    #' Removes from the EEG data the samples of any epoch in `epochs`.
    #'
    #' @param epochs (`int`)
    #' @return void
    drop_epochs = function(epochs) {
        df <- self$data 
        self$data <- droplevels(self$data[!self$data$Epoch %in% epochs, ])[-2]
    },

    #' @description
    #' Removes from the EEG data the samples of any epoch-subepoch pair
    #' in `epochs`, `subepochs`.
    #' @param epochs (`int`)
    #' @param subepochs (`int`)
    #' @return void
    drop_subepochs = function(epochs, subepochs) {
        contaminated <- as.factor(paste(epochs, subepochs))
        df <- self$data
        df <- df %>% tibble::add_column(
            Pairs = as.factor(paste(df$Epoch, df$Subepoch)),
            .after = "Subepoch"
        )
        df <- droplevels(df[!df$Pairs %in% contaminated, ])[-c(4)]
        self$data <- df
    },

    # ---- Artifact related functions

    #' @description
    #' Performs MV-CAPA analysis on the EEG data; formats the results
    #` and sets them to the `canoms` and `panoms` fields.
    #'
    #' @param alpha (`numeric`)\cr 
    #' Only collective anomalies with strength greater than `alpha` are kept.
    #' Defaults to `8`.
    #' @param beta (`numeric`)\cr 
    #' Only point anomalies with strength greater than `beta` are kept.
    #' Defaults to `1`.
    #' @return void
    artf = function(alpha = 0.1) {
        print("Starting artifact analysis. This may take a few minutes...")
        analysis <- anomaly::capa(self$data[, -c(1:3)], type = "mean")
        canoms <- anomaly::collective_anomalies(analysis)
            tibble::as_tibble()
        panoms <- anomaly::point_anomalies(analysis) %>%
            tibble::as_tibble()
        self$canoms <- canoms
        self$panoms <- panoms
        self$sfilter(alpha)
        self$set_anom_time_features()
    },

    #' @description
    #' Performs MV-CAPA analysis on the EEG data by steps of size `step_size`
    #' (in seconds); formats the results and sets them to the `canoms` and
    #' `panoms` fields.
    #'
    #' @param step_size (`int`) Seconds per analyzed step.
    #' @param alpha (`numeric`) The result will exclude anoms with normalized strength <= alpha.
    #' The value of `alpha` must satisfie 0 ≤ alpha ≤ 1.
    #' Defaults to `0.05`.
    #' @param type (`string`) Should be "mean" or "meanvar", depending on 
    #' whether to test for changes in mean or in mean and variance.
    #' @param verbose (`bool`) Print log messages of the analysis process?
    #' 
    #' @return void
    artf_stepwise = function(step_size = 30, alpha = 0, type = "mean",
                             verbose = FALSE) {
        print("Starting analysis. This may take a couple of minutes...")
        exclude <- c("Time", "Epoch", "Subepoch")
        q <- self$data$Time %/% step_size
        t <- self$data %>%
            tibble::add_column(AnRegion = as.factor(q + 1), .after = 1)
        grouped <- dplyr::group_by(t[, !names(t) %in% exclude], AnRegion) %>%
            dplyr::group_map(~ anomaly::capa(x = scale(.x), type = type))
        gc()
        canoms <- grouped %>% lapply(function(x) {
                                         anomaly::collective_anomalies(x)
                                    })
        gc()
        panoms <- grouped %>% lapply(function(x) anomaly::point_anomalies(x))
        gc()
        mps <- as.vector(table(t$AnRegion)) # Measures per step.
        self$set_anom_dfs(mps, panoms, canoms, alpha)
    },

    #' @description
    #' A helper function to be privately used in artf_stepwise. Transforms each 
    # ` of two lists of data frames, `panoms` and `canoms`, into a unique 
    #' tibble with the joined data of each analysis; it adjusts the location
    #` features of these data frames to match the full `$data` instead of 
    #` simply the step to which they correspond.
    #'
    #' @param mps (`vector`) A vector containing the number of measures per analyzed step of data.
    #' @param panoms (`list`) List of data frames, each comming from calling
    #' anomaly::point_anomalies on an analyzed subset of the EEG data.
    #' @param canoms (`numeric`) List of data frames, each comming from calling
    #' anomaly::collective_anomalies on an analyzed subset of the EEG data.
    #' @param alpha See artf_stepwise which transfers its alpha parameter to this.
    #' @return void
    set_anom_dfs = function(mps, panoms, canoms, alpha=0.05) {
        mps <- c(0, mps[-length(mps)])
        self$canoms <- mapply(function(x, y) {
                            x %>% dplyr::mutate(
                                                start = start + sum(mps[1:y]),
                                                end = end + sum(mps[1:y]))
                            },
                            canoms,
                            seq_along(canoms), SIMPLIFY = FALSE) %>%
            bind_rows() %>%
            as_tibble()

        self$panoms <- mapply(function(x, y) {
                         x %>% mutate(location = location + sum(mps[1:y]))
                        },
                        panoms, seq_along(panoms),
                        SIMPLIFY = FALSE) %>%
            bind_rows() %>%
            as_tibble()
        self$set_anom_time_features()
        self$sfilter(alpha)
    },

    #' @description
    #' Private wrapper.
    #' @return void
    set_anom_time_features = function() {
        self$panoms$Time <- self$data$Time[unlist(self$panoms[1])]
        self$canoms$Time <- self$data$Time[unlist(self$canoms[1])]
        self$canoms <- set_epochs(self$canoms, subepochs = TRUE)
        self$panoms <- set_epochs(self$panoms, subepochs = TRUE)
    },

    #' @description
    #' Private wrapper.
    #' @return void
    get_contaminated_channels = function() {
        chans <- union(
            unique(self$canoms$variate),
            unique(self$panoms$variate)
        )
        return(chans)
    },

    #' @description
    #' (Private)
    #' Returns a df whose format is idoneous for plotting the analysis results.
    #' Such format consists of the EEG Time measure (x axis), the raw
    #' EEG data of the chanel and a column that is filled with NaN except
    #' at each time instance where any type of anomaly was found.
    #'
    #' @param chan (`int`) Integer in the range [1, number_of_eeg_signals].
    #' @param s Plot should begin at epoch s.
    #' @param e Plot should end at epoch e.
    #' @return data.frame
    set_plot_data = function(chan, s = 0, e = 0) {
        data <- self$data
        canoms <- dplyr::filter(self$canoms, variate == chan)
        panoms <- dplyr::filter(self$panoms, variate == chan)
        subsetting <- FALSE
        if (s != 0 && e != 0) {
            if (!(s %in% data$Epoch) || !(e %in% data$Epoch)) {
                stop("Invalid epoch bounds")
            }
            data$ExplicitIndex <- seq_len(nrow(data))
            canoms <- dplyr::filter(canoms, Epoch %in% c(s:e))
            panoms <- dplyr::filter(panoms, Epoch %in% c(s:e))
            data <- dplyr::filter(data, Epoch %in% c(s:e))
            subsetting <- TRUE
        }
        # Get all indexes between start and end of canoms
        locations <- mapply(function(x, y) x:y, canoms$start, canoms$end)
        # Unite with point anomalies
        locations <- union(unlist(locations), unlist(panoms$location)) %>%
                    as.integer()
        # When working with the original 'self$eeg$data', 'locations' matches
        # anomalous sample indexes. If dealing with a subset, 'locations'
        # corresponds to 'self$eeg$data' indexes, not the subset.
        # These cases are handled separately. While using the 'ExplicitIndex'
        # column for both cases is an option, it's more efficient to check a
        # boolean condition than to directly add a potentially large column.
        if (!subsetting) {
            time_of_anomalies <- unlist(data[locations, 1]) %>%
                                lubridate::as_datetime()
            values <- unlist(data[locations, chan + 3])
        }else {
            anom_samples <- dplyr::filter(data, ExplicitIndex %in% locations)
            time_of_anomalies <- unlist(anom_samples[, 1]) %>%
                                lubridate::as_datetime()
            values <- unlist(anom_samples[, chan + 3])
        }
        df <- tibble::tibble(X = time_of_anomalies, Y = values)
        return(df)
    },

    #' @description
    #' Plots analysis results for a specific channel.
    #'
    #' @param chan An integer that points to the EEG channel.
    #' @param s Plot should begin at epoch s.
    #' @param e Plot should end at epoch e.
    #' @param size (optional) The size of the red dots that signal an anomaly.
    #' @return ggplot
    plot_channel_artifacts = function(chan, s = 0, e = 0, size = 0.2) {
        df <- self$set_plot_data(chan, s, e)
        eeg <- self$plot_channel(chan, s, e)
        p <- eeg +
            ggplot2::geom_point(
                data = df, ggplot2::aes(X, Y),
                inherit.aes = FALSE, color = "red",
                size = size
            )
        return(p)
    },

    #' @description
    #' Plots all channels with any kind of anomalies from epoch `s` to `e`.
    #' Mark anomalies with dots of radius `size`.
    #' @param s Starting epoch
    #' @param e Ending epoch
    #' @param size Size of the red dots that mark anomalies in the plot
    #' @return cowplot
    plot_artifacts = function(s = 0, e = 0, size = 0.2) {
        plots <- list()
        channels <- self$get_contaminated_channels()
        for (channel in channels) {
            p <- self$plot_channel_artifacts(channel, s, e, size)
            plots[[channel]] <- p
        }
        return(cowplot::plot_grid(plotlist = plots, align = "v", ncol = 1))
    },

    #' @description
    #' Saves to `dir` a plot of all artifacts. Each plot contains a number of
    #' of epochs determined by `epochs_per_plot`. If `ignore_clean` is `TRUE`,
    #' then any segment not containing anomalies is not plotted.
    #' 
    #' @param epochs_per_plot Size of the x-axis in epochs.
    #' @param dir Directory where plots will be saved.
    #' @param ignore_clean Skip plotting clean segments?
    #' @return void
    gen_artifact_plots = function(epochs_per_plot,
                                  dir = getwd(),
                                  ignore_clean=TRUE) {
        epochs_per_plot <- epochs_per_plot - 1
        s <- as.numeric(head(self$data$Epoch, 1))
        end <- as.numeric(tail(self$data$Epoch, 1))
        while (s + epochs_per_plot <= end) {
            if (ignore_clean && !(self$has_anoms(s, s + epochs_per_plot))) {
                s <- s + epochs_per_plot + 1
                next
            }
            p <- self$plot_artifacts(s, s + epochs_per_plot)
            fname <- paste(dir, "/", s, ".png", sep = "")
            ggplot2::ggsave(fname, p)
            s <- s + epochs_per_plot + 1
        }
    },

    #' @description
    #' Determines whether a data segment contains anomalies. If `e` parameter 
    #' is zero (default), the function test whether any anomaly exists in 
    #' the whole EEG data.
    #'
    #' @param s Starting epoch of the segment
    #' @param e Ending epoch of the segment, defaults to zero.
    #' @return bool
    has_anoms = function(s, e = 0) {
        if (e == 0) {
            e <- s
        }
        canoms <- dplyr::filter(self$canoms, Epoch %in% c(s:e))
        panoms <- dplyr::filter(self$panoms, Epoch %in% c(s:e))
        l <- list(canoms, panoms)
        i <- which(lapply(l, nrow) == 0)
        # If both data frames have length 0
        if (length(i) > 1) {
            return(FALSE)
        }
        return(TRUE)
    },


    #' @description
    #' Normalizes the strengths of point and collective anomalies
    #' given a normalizing function `f`, which defaults to min-max
    #' normalization, and removes anomalies whose strength is inferior
    #' to `x`.
    #'
    #' @param x Threshold value
    #' @param f Normalization function (defaults to min-max).
    #' @return void
    sfilter = function(x, f = minmax_normalization) {
        self$canoms$mean.change <- f(self$canoms$mean.change)
        self$panoms$strength <- f(self$panoms$strength)
        self$canoms <- self$canoms %>% dplyr::filter(mean.change >= x)
        self$panoms <- self$panoms %>% dplyr::filter(strength >= x)
    },

    #' @description
    #' COMPLETE
    #'
    #' @return void
    get_contaminated_epochs = function() {
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

    #' @description
    #' Performs artifact rejection by means of dropping from the EEG data
    #' all segments corresponding to epoch-subepoch pairs
    #' that were found to contain anomalies.
    #'
    #' @return void
    artf_reject = function() {
        epoch_data <- self$get_contaminated_epochs()
        clone <- self$clone()
        clone$drop_subepochs(epoch_data$Epoch, epoch_data$Subepoch)
        clone$canoms <- clone$panoms <- tibble::tibble()
        return(clone)
    },

    # --- PSD ---


    #' @description
    #' Plots the spectrogram of an EEG signal.
    #' @param channel Integer pointing to an EEG signal.
    #' @param max_freq Maximum frequency of the plot.
    #' @param freq Aggregate frequency to set plot resolution.
    #'
    #' @return void
    spectrogram = function(channel, max_freq = 30, freq = 4) {
        p <- rsleep::spectrogram(unlist(self$data[, -c(1:3)][channel]),
                            sRate = self$fs, maxFreq = max_freq,
                            freq = freq)
        return(p)
    },

    #' @description
    #' Computes the PSD of an EEG signal. It is returned in data frame format.
    #' @param channel Integer pointing to an EEG signal.
    #'
    #' @return data.frame
    channel_psd = function(channel) {
        welch <- gsignal::pwelch(as.matrix(self$data[, -c(1:3)][channel]), fs = self$fs)
        psd <- welch$spec %>%
            apply(log10, MARGIN = 2) %>%
            tibble::as_tibble()
        psd$Fqc <- welch$freq
        return(psd)
    },

    #' @description
    #' Computes the PSD of each EEG signal. Result is set to the $psd field.
    #'
    #' @return data.frame
    compute_psd = function() {
        pwelch <- gsignal::pwelch(as.matrix(self$data[, -c(1:3)]), fs = self$fs)
        psd <- pwelch$spec %>%
            apply(log10, MARGIN = 2) %>%
            tibble::as_tibble()
        psd$Fqc <- pwelch$freq
        self$psd <- psd
    },

    #' @description
    #' Produces a ggplot of the PSD. `compute_psd` must have been called before.
    #'
    #' @param xlim Maximum frequency to show
    #' @return ggplot
    plot_psd = function(xlim = 250) {
        tall <- reshape2::melt(self$psd, id.vars = "Fqc")
        p <- ggplot2::ggplot(tall, ggplot2::aes(Fqc, value, col=variable)) +
            ggplot2::geom_line() +
            ggplot2::xlim(c(0, xlim))
        return(p)
    },

    #' @description
    #' Produces a plotly plot of the PSD. `compute_psd` must have been called before.
    #'
    #' @param xlim Maximum frequency to show
    #' @return plotly figure
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
    },

    #' @description
    #' Performs zero-mean centering per epoch. Useful pre-processing step for
    #' certain analyses (e.g. ICA).
    #' @return void
    center_channels = function() {
        x <- eeg$data[, -c(1, 3)]
        centered <- x %>%
                    tibble:group_by(Epoch) %>%
                    tibble::group_modify(~tibble::as_tibble(scale(.x))) %>%
                    tibble::as_tibble(.name_repair = "minimal")
     self$data[, -c(1:3)] <- centered[, -c(1)] #1st col of centered = Epoch
    },

    #' @description
    #' Performs spindle detection either on a specific signal (if `channel` is
    #' non-zero) or for each EEG signal. The two available methods are
    #' Sigma Index and Relative Spindle Power. For details on these, consult the
    #' GitHub README. The result is set to the $spindles field.
    #' '
    #'
    #' @param channel On which signal to detect spindles? If zero (default),
    #'                all signals are subjected to spindle detection.
    #' @param method Either "sigma_index" (default) or "rsp".
    #' @param filter A boolean determining whether to filter spindles according
    #'                to standard strength thresholds.
    #' @return void
    spindle_detection = function(channel = 0, # channel = 0 -> whole EEG
                                 method = "sigma_index",
                                 filter = TRUE) {
        if (method == "sigma_index") {
            f <- sigma_index
            threshold <- 4.5
        }else if (method == "rsp") {
            f <- relative_spindle_power
            threshold <- 0.22
        }else {
            stop("Invalid spindle detection method. `method` argument should be:
            'rsp': for Relative Spindle Power detection, or
            'sigma_index': for Sigma Index detection.")
        }
        if (channel != 0) {
            data <- self$data[c(1, channel + 3)] %>% window_eeg_data(1)
        }else {
            data <- self$data[, -c(2:3)] %>% window_eeg_data(1)
        }
        result <- data %>%
            dplyr::group_by(Group = Windows) %>%
            dplyr::summarize(dplyr::across(-c(Time, Windows), \(x) f(x, fs = self$fs))) %>%
            dplyr::mutate(Group = as.numeric(stringr::str_extract(Group, "\\d+\\.?\\d*?(?=,)"))) %>%
            dplyr::rename(Second = Group) %>%
            na.omit()
        if (filter) {
            result <- result %>%
                filter(across(-Second, ~ . > threshold)
                %>% rowSums() > 0) %>%
                dplyr::mutate(across(-Second, ~ ifelse(. < threshold, 0, .))) %>%
                dplyr::select(where(~ any(. != 0)))
        }
        self$spindles <- result
    },


    #' @description
    #' Plots the distribution of spindles across the EEG. `spindle_detection`
    #' must have been previously called.
    #'
    #' @param channel On which signal to detect spindles? If zero (default),
    #'                a cumulative index of the spindles on all signals is
    #'                plotted.
    #' @param time_axis Show distribution over "epoch", "second", "minute" or "hour"?
    #' @param xbins Size of the x-bins in the plot.
    #' @param ybins Size of the y-bins in the plot.
    #' @param from Start of the plot. Defaults to 0 (whole EEG).
    #'
    #' @return void
    plot_spindle_distribution = function(channel = 0,
                                         time_axis = "epoch",
                                        xbins = 10,
                                        ybins = 10,
                                        from = 0) {
        time_resolution <- c(second = 1 / 60,
                            minute = 1 / (60 ^ 2),
                            hour = 1 / (60 ^ 3),
                            epoch = 1 / 30)
        data <- self$spindles
        if (channel == 0) {
            ylabel <- "Cumulative Index Score"
            y <- rowSums(data[, !colnames(data) %in% c("Second")])
            plot_title <- "Cumulative spindle distribution over time"
        }else{
            y <- data[[channel + 1]]
            ylabel <- "Index Score"
            plot_title <- paste(colnames(data)[channel + 1],
                                ": Spindle distribution over time",
                                sep = "")
        }

        fig <- plotly::plot_ly(self$spindles,
            x = ~(Second * time_resolution[time_axis]),
            y = y)
        fig <- fig %>%
          plotly::add_histogram2d(
                nbinsx = xbins,
                nbinsy = ybins,
                ybins = list(start = from)
            ) %>%
          plotly::layout(
            title = plot_title,
            xaxis = list(title = paste(time_axis, "s", sep = "")),
            yaxis = list(title = ylabel)
          )
        return(fig)
    }
    )
)
