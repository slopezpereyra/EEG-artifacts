
# EEG Automated Artifact Detection

  
  
:last_quarter_moon_with_face:
> :microscope: This software is an implementation of [Multivariate-CAPA](https://arxiv.org/abs/1806.01947) for EEG artifact detection.

  

- [EEG Automated Artifact Detection](#eeg-automated-artifact-detection)

- [Overview](#overview)

- [Installation](#installation)

- [Loading EEG data](#loading-eeg-data)

- [Handling EEG data](#handling-eeg-data)

- [Anomaly detection](#anomaly-detection)

- [Example analyses](#example-analyses)

- [1. Eleventh minute](#1-eleventh-minute)

- [B. Seventeenth minute](#b-seventeenth-minute)

- [Analyzing the whole EEG](#analyzing-the-whole-eeg)

  
### Installation

  

In order to perform automated artifact detection, we need to install the `artifactor` package in R. Simply run the following commands in the R console.

  

```r
library(devtools)
install_github("slopezpereyra/EEG-artifacts")
library(artifactor)
```
  

### Loading EEG data

  

To load EEG data, we use the `load_eeg(data_file, signals_file = NULL)` function. This function takes a data file and an optional signals file as arguments.

  

> :bulb: Exporting EDF files via _EDF Browser_ in a .csv format will always yield a data file and a signals file (among others). The optional signals file contains information pertaining to the EEG channels, and its inclusion allows the package to set appropriate channel names.

  

For example, after previously exporting an EDF file named `test` in a .csv format,

  

```r
eeg <- load_eeg("test_data.txt",
				"test_signals.txt") %>%
				na.omit()
View(eeg@data)
```

  

![EEG data](https://i.ibb.co/M9CqQzG/Screenshot-from-2022-09-04-16-20-50.png)

  

The `load_eeg` function returns an `eeg` object containing the `@data` and `@signals` atrributes, the latter being an empty data frame if no signals file was provided.

We can consult the sampling frequency if we ignore it.

```
get_sampling_frequency(eeg)
> 500 # In Hz
```



  

### Handling EEG data

  

Due to the size of EEG data, it may be desirable to work with only a subset of the record, or a low resolution version of it. This is accomplished via the `subset_eeg(eeg, start, end)` and `resample(eeg, n)` functions.

`subset_eeg` is a sub-setting function that extracts a time period of data in seconds. `resample` reduces the EEG data by keeping only one very $n$ observations.
 
For example,
  

```r
tenth_min_data <  subset_eeg(eeg, 60  *  10, 60  *  11)
resampled_eeg < resample(tenth_min_data, 20)
```

  

defines two new `eeg` objects. The first contains only the tenth minute of the original record; the latter contains that same timespan with a resolution such that $5$ entries (instead of $500$) make up a second. If we plot the F3-A2 channel of both sets, the difference is evident.

  

![Full resolution](https://i.ibb.co/PgP1S3P/plot.png)

  

Artifact rejection analysis may be conducted on any of the two formats, with an expectable accuracy vs. performance trade-off.

  

Keep in mind that lowering the data resolution is not mathematically equivalent to performing a low-pass filter. If one should want to apply a low-pass filter, then the `low.pass(eeg, n)` function should be called, where a $\frac{1}{n} \text{Hz}$ [Butterworth filter](https://en.wikipedia.org/wiki/Butterworth_filter) is applied to all channels. For example,

```r
minute <- subset_eeg(eeg, 60, 120) # Second minute of record
fminute <- low_pass(fminute, 20) # Apply 20Hz filter
```
 
defines two `eeg` objects containing the second minute of record, the latter with a low-pass filter. Visually, they differ as shown in the image below.

![enter image description here](https://i.ibb.co/HnG5jTc/plot-3.png)

  

> :bulb: Note that, while performing artifact detection on low-resolution data is faster, that is not the case with low-pass filtered data. The reason is that `lower.res(eeg, n)` reduces the amount of values in the record by a factor of $\frac{1}{n}$, while the `low.pass` returns an `eeg` object with the same number of records. In other words, `lower.res` is a subsetting function, while `low.pass` is a matrix transformation. The degree to which low-resolution scaling negatively impacts artifact detection is still to be determined.

  

### Anomaly detection

Anomaly detection is carried out via the CAPA statistical method ([Fisch, Eckley & Fearnhead, 2021](https://onlinelibrary.wiley.com/doi/full/10.1002/sam.11586)). In order to adapt CAPA to the specificities of EEG data, the package conducts a substantial amount of operations. However, these are carried out _sub rosa_, so that performing an analysis is straightforward. Simply call the `analyze(df, start, end, ...)`. :microscope:

  

The `analyze` function returns an `analysis` object. An `analysis` object has

  

- a `@canoms` attribute: data frame of all collective anomalies;

- a `@panoms` attribute: data frame of all point anomalies;

- an `@eeg` attribute: containing the eeg upon which the analysis was conducted.

  

The `analysis` class has numerous plotting methods, the most important of which is `plot(analysis)` (used to create the graphs of the forthcoming section).

  

#### Example analyses

  

We will now choose certain portions of our test EEG record that seem to contain more or less subtle artifacts and perform an analysis on them.

  

##### 1. Eleventh minute

  

Firstly, we will perform anomaly detection on the eleventh minute of data with $\alpha = 8$.

  

> :bulb: Alpha is the threshold value determining how far from the estimated distribution mean a sequence ought to stray to be considered anomalous. Experimentation has found a good standard value is $8$.

  

```r

analysis <- analyze(eeg, 60  *  11, 60  *  12, alpha  =  8)
plot_analysis(analysis)

```

  

![Analysis](https://i.ibb.co/7KgzB77/analysis.png)

  

The algorithm is finding a long artifact in the last ten seconds of data, while all channels also show anomalies around minute 11:15. If we are curious to see a particular channel, we can call a specific plot for it. Let us see what is going on with LOC-A2, the seventh channel.

  

```r
plot_anoms(analysis, channel  =  7)
```

  

![LOC-A2](https://i.ibb.co/DgrQH7G/analyisis-c7.png)

  

##### B. Seventeenth minute

  

```r
analysis <- analyze(eeg, 17  *  60, 18  *  60, alpha  =  8)
plot(analysis)
```

  

![enter image description here](https://i.ibb.co/zQDt61k/plot-4.png)

  

## Analyzing the whole EEG

  

Our package uses the `setpwise_analysis(eeg, step_size, res, ...)` function to conduct stepwise artifact analysis on large EEGs.

  

> :bulb: By stepwise analysis we mean the separate analysis of the sequences $a_0, a_1, a_2, ..., a_n$, where $a_i$ contains the values of the EEG from seconds $i \cdot s$ to $i(s+1)$, with $s =$ `step_size`.  For an explanation of why stepwise (instead of direct) analysis is required, and the logic behind the algorithm, see [this document](https://github.com/slopezpereyra/EEG-artifacts/blob/main/docs/stepwise.md).

 
For example, `setpwise_analysis(eeg, 30, res = 1)` would perform artifact detection over the raw EEG by analyzing each of its thirty second epochs separately.

The function returns a single analysis object containing the joint anomalies found across the whole analyzed record. 

We will exemplify the function with a small subset of the data; namely, the first five minutes of record of our test EEG. During this period, the subject was requested to perform certain motor activities, like gnashing his teeths and moving his legs. Because of this, the record contains a good deal of anomalous segments. These anomalies, although not artifacts —for they were purposeful—, serve us good for showcasing purposes.

```r
first_five_mins <- subset_eeg(eeg, 0, 5  *  60)
analysis <- setpwise_analysis(first_five_mins, 60, res  =  1)
```
The console outputs the time it took for the analysis to complete. In this case, the five minutes were done in $\approx10.5$ seconds., although this could have been faster had we scaled the resolution. 

We can plot the fifth minute of this analysis, for example, just as we plotted analysis objects before.

![enter image description here](https://i.ibb.co/t8wRyLp/0-4-0-to-0-5-0.png)

