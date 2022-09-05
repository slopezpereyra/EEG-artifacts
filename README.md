﻿# EEG Automated Artifact Detection

> :microscope: This software is an implementation of [Multivariate-CAPA](https://arxiv.org/abs/1806.01947) for EEG artifact detection.

### Overview

The aim of this document is to succinctly show how our automated artifact detection package works. This should be sufficient to allow other researchers to begin using the package.

### Installation

In order to perform automated artifact detection, we need to install our `artifactor` package in R. Simply run the following commands in the R console.

    library(devtools)
    install_github("slopezpereyra/EEG-artifacts")
    library(artifactor)

And that is it! Our package is ready to go. :fire::smiling_imp:

### Loading EEG data

To load EEG data, we use the `load_eeg(data_file, signals_file)` function and pass data and signal files as arguments.

> :bulb: Exporting EDF files via _EDFBrowser_ in a .csv format will always yield data and a signal files. This are the files required by the function. The signals file contains information pertaining to the EEG channels (e.g., their labels and their min and max values), while the data file holds the actual EEG measures.

For example, after previously exporting an EDF file named `test` in a .csv format,

    data <- load.eeg("test_data.txt",
    				"test_signals.txt") %>%
    		na.omit()
    View(data)

![EEG data](https://i.ibb.co/M9CqQzG/Screenshot-from-2022-09-04-16-20-50.png)
Channel labels are automatically recognized and appropriately set as column names.

### Handling EEG data

EEG data is extremely large. For that reason, it is often desirable to work with only a subset of the record, or a low resolution version of it. This is accomplished via the `partition.eeg(df, start, end)` and `lower.res(df, n)` functions. The latter is equivalent to applying a low-pass filter, since it reduces the number of events per unit of time (only one every $n$th values are kept). The first is a sub-setting function that extracts a time period of data in seconds.

For example,

```js
tenth_min_data < -partition.eeg(data, 60 * 10, 60 * 11);
lowpass < -lower.res(tenth_min_data, 20);
```

defines two new data sets. The first contains the raw tenth minute of the record; the latter, that same timespan with a resolution such that $5$ entries (instead of a hundred) make up a second. If we plot the F3-A2 channel of both sets, the difference is evident.

![Full resolution](https://i.ibb.co/PgP1S3P/plot.png)

We could conduct artifact rejection analysis on any of the two formats, with the expectable accuracy vs. performance trade-off.

### Anomaly detection

Anomaly detection is done via the CAPA statistical method ([Fisch, Eckley & Fearnhead, 2021](https://onlinelibrary.wiley.com/doi/full/10.1002/sam.11586)). In order to adapt CAPA to the specificities of EEG data, the package conducts a substantial amount of operations. However, this are carried out _sub rosa_, so that performing an analysis is straightforward. Simply call the `analyze(df, start, end, ...)`. :microscope:

The `analyze` function returns an `analysis` object defined by a data frame of all collective anomalies (or anomalous segments of data), a data frame of all point anomalies, and the original data where the analysis was performed. Results can be directly plotted via suittable functions.

#### Example analysis

We will perform anomaly detection on the tenth minute of data with $\alpha = 8$.

> :bulb: Alpha is the threshold value determining how far from the estimated distribution mean a sequence ought to stray to be considered anomalous. Experimentation has found a good standard value is $8$.

    analysis <- analyze(data, 60 * 11, 60 * 12, alpha = 8)
    plot(analysis)

![Analysis](https://i.ibb.co/7KgzB77/analysis.png)

The algorithm is finding a long artifact in the last ten seconds of data, while all channels also show anomalies around minute 11:15. If we are curious to see a particular channel, we can call a specific plot for it. Let us see what is going on with LOC-A2, the seventh channel.

    draw.anomalies(analysis, channel = 7)

![LOC-A2](https://i.ibb.co/DgrQH7G/analyisis-c7.png)

## Analyzing the whole EEG

Lastly, our package provides the `analyze.stepwise(df, step_size, res, ...)` function. This function performs stepwise analysis on the whole EEG data frame passed as argument.

> :bulb: By stepwise analysis we mean the separate analysis of the sequences $a_0, a_1, a_2, ..., a_n$, where $a_i$ contains the values of the EEG from seconds $i \cdot s$ to $i(s+1)$, with $s =$ `step_size`.

For example, `analyze.stepwise(data, 30, res = 1)` would perform artifact detection over the raw EEG by analyzing each of its thirty second epochs separately.

The `analyize.stepwise` function saves the analysis plots of each separate analysis as `.png` files and also saves a single `.csv` spreadsheet containing each epoch-subepoch pair containing anomalies.

Because the output of this function is large (dozens of `,png` images images and a `.csv` file), we will not show it here. But this is the function one should use when aiming at analyzing not a specific portion of the EEG, but its entirety.
