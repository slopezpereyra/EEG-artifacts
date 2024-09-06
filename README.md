
# EEG Computational Toolkit
  

> :last_quarter_moon_with_face: Developed at the [Laboratory for the Study of
> Sleep Slow-wave activity](https://www.med.upenn.edu/slowwavelab/)

A scientific package for computational EEG analysis.

- [EEG Computational Toolkit](#eeg-computational-toolkit)
    - [Installation and import](#installation-and-import)
    - [Loading EEG data](#loading-eeg-data)
    - [The `EEG` object](#the-`eeg`-object)
    - [EEG Visualization](#eeg-visualization)
    - [Resampling and resolution scaling](#resampling-and-resolution-scaling)
    - [Filtering](#filtering)
    - [Artifact detection](#artifact-detection)
  - [Power spectrum analysis](#power-spectrum-analysis)
  - [Example artifact detection](#example-artifact-detection)

This package has a [Julia
alternative](https://slopezpereyra.github.io/EEGToolkit.jl/dev/) and has been
deprecated in favor of it.

### Installation and import

Install the `remotes` package with `install.packages('remotes')` and run in R

```R 
remotes::install_github("slopezpereyra/EEG-toolkit")
```

*Note:* If some of the packages on which this package depends (e.g.
`tidyverse`) cannot be installed due to unmet dependencies, try running 

```
apt-get update
# System dependencies for ALL dependencies
sudo apt-get install -y \
    libgdal-dev gdal-bin libproj-dev proj-data proj-bin libgeos-dev \
    libcurl4-openssl-dev libssl-dev libfontconfig1-dev libxml2-dev  \ 
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

and then repeat the `install_github` command. Most of these libraries are 
probably installed in your system already.

### Loading EEG data

To load EEG data, we use the `EEG$new(data_file, ...)` initiation function. 
This function takes the path (`string`) to an EDF or CSV file as argument and
initializes an `EEG` object with the data read. 

Reading an EDF takes more or less the same time than reading a CSV, but the
resulting data is more polished. For example, because EDF files contain channel
names, column names will be immediately set when reading this format, while they
must be manually set if reading .csv files.


For example, if we have an `eeg.csv` file in some path relative to our working
directory, we call

```r 
eeg <- EEG$new("/relative/path/eeg.csv") # could be eeg.edf as well...
```

### The `EEG` object

The `EEG$new` function returns an `EEG` R6 object with the following fields: 

- `$data` (`tibble`): The EEG data as a data frame.
- `$canoms` (`tibble`): A data frame with collective anomalies data. If artifact
  detection was not yet performed, it defaults to an empty data frame.
- `$panoms` (`tibble`): A data frame with point anomalies data. If artifact
  detection was not yet performed, it defaults to an empty data frame.
- `$psd` (`tibble`) : The estimated power spectrum density of the EEG. If PSD computation
  was not yet performed, it defaults to an empty data frame.
- `$spindles` (`tibble`): The spindles detected by any of the automated spindle detection
  algorithms, it defaults to an empty data frame.
- `$fs` (`int`): The sampling frequency $f_s$ of the EEG.

Most methods of the EEG class work in-place. For example, take the
following lines. 

``` r
eeg <- EEG$new("/relative/path/eeg.csv") 
eeg$subset(1, 10)  
```

This subsets the `eeg` object from epoch $1$ to $10$ *in-place*, meaning that the
original EEG was modified. For that reason, it is advisable to create a copy of
the original EEG object before proceeding with any analysis: 

``` r
eeg <- EEG$new("/relative/path/eeg.csv")
safe_clone <- eeg$clone() 
eeg$subset(1, 10)
```

Now, albeit we are modifying the `eeg` object, we always have our original EEG
in the `safe_clone` object.

### EEG Visualization

We can produce static or interactive plots to get a visual sense of our record.
For example,

```r 
eeg$subset_by_seconds(900, 960) # Subset by time, not epochs
eeg$plot()
``` 

![enter image description here](https://i.ibb.co/0X4GG8T/Screenshot-from-2022-12-05-13-12-17.png)

We could have also produced an interactive plot to inspect our EEG record, or a
portion of it, live. For example,

```r 
eeg$iplot()
```

produces the following plot.


![Alt Text](https://i.ibb.co/0XgxhKv/ezgif-1-2c5fd0d1e6.gif)

  

### Resampling and resolution scaling

EEG data is large. For a sampling rate $f_s = 500$, a single $30$ seconds epoch
contains $15.000$ observations.  Hence, it is often desirable to work either
with subsets of the record, or lower resolution version of it. These package
provides functions to make subsetting and resolution-scaling practical and easy.

Subsetting methods have already been showcased above. Resampling, on the other
hand, is performed  via the `resample(n)`. Resampling will keep only one
every $n$ observations, producing a lower resolution version of the record. This
can help accelerate different types of analysis, such as artifact detection and
power spectrum analysis, as well as static and interactive plotting.
 
### Filtering 

Low-pass, high-pass and bandpass filters are available via the `low_pass,
high_pass` and `bandpass` functions. For example, 


```r 
eeg <- EEG$new("some/relative/path/eeg.csv")
eeg$subset_by_seconds(60, 120) # Reduce the EEG to the second minute of record 
eeg$plot()
eeg$low_pass(20) # Apply 20Hz filter
eeg$plot() # Plot once more.
```
 
We have plotted the EEG before and after applying a low-pass filter. The
resulting plots differ as shown below.

![enter image description here](https://i.ibb.co/HnG5jTc/plot-3.png)


### Artifact detection

Artifact detection is carried out via the CAPA statistical method ([Fisch,
Eckley & Fearnhead,
2021](https://onlinelibrary.wiley.com/doi/full/10.1002/sam.11586)). CAPA is
adapted to the specificities of EEG data via _sub rosa_ operations. To perform
artifact detection simply call `artf(eeg, start, end, ...)`. :microscope:

Once an analysis is performed, the `$canoms` and `$panoms` attributes of the EEG
are filled with the analysis results. The first pertains to collective
anomalies, the second to point anomalies. 

Plotting functions allow us to observe the results after artifact detection was
performed. For example.

``` r
eeg <- EEG$new('some/relative/path/eeg.csv')
eeg$subset_by_seconds(900, 960)
eeg$artf() # Performs direct artifact detection
eeg$plot_artifacts()
```

![enter image description
here](https://i.ibb.co/BjL6fDR/Screenshot-from-2022-12-05-13-18-15.png)

Marked in red are the primary suspects.  The fourth channel is omitted, since no
anomaly was found. 

It should be noted that anomalies have different strengths. We can easily filter
an analysis object so as to keep only anomalies stronger than a certain
threshold value. For example,

``` r
eeg$sfilter(0.1)    # After min-max normalization, 
                    # keep only anomalies with strength >= 0.1 
eeg$plot_artifacts()
```

![enter image description
here](https://i.ibb.co/djz0v74/Screenshot-from-2022-12-05-13-21-10.png)

The complexity of the algorithm depends primarily on the complexity of `capa`,
which is detailed in Fisch, Eckley and Fearnhead. On a 15.5 million samples EEG
record with 8 channels, corresponding to a full night of sleep, artifact
detection with a `step_size` parameter of `30 * 5` seconds (each five epoch
period was independently analyzed ) took on average ≈ 8.6 minutes. The
specifications for one of these tests is given below for reference.

```
========================================================================= 
R Version: 4.2.2 Patched (2022-11-10 r83330) 
Operating System: Linux 6.2.0-34-generic #34-Ubuntu SMP PREEMPT_DYNAMIC Mon Sep 4 13:06:55 UTC 2023 
Base Packages: stats graphics grDevices utils datasets methods base Other Packages: forcats_1.0.0 stringr_1.5.0 purrr_1.0.1 readr_2.1.4 tidyr_1.3.0 tibble_3.2.1 ggplot2_3.4.2 tidyverse_1.3.2 dplyr_1.1.2 logr_1.3.4 
========================================================================= 

# A tibble: 15,562,500 × 11
    Time Epoch Subepoch `F3-A2` `F4-A1` `C3-A2` `C4-A1` `O1-A2` `O2-A1` `LOC-A2`
   <dbl> <fct> <fct>      <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>
 1 0     1     1        0.00478 0.00478 0.00478 0.00478 0.00478 0.00478  0.00478
 2 0.002 1     1        0.00478 0.00478 0.00478 0.00478 0.00478 0.00478  0.00478
 3 0.004 1     1        0.00478 0.00478 0.00478 0.00478 0.00478 0.00478  0.00478
 4 0.006 1     1        0.00478 0.00478 0.00478 0.00478 0.00478 0.00478  0.00478
 5 0.008 1     1        0.00478 0.00478 0.00478 0.00478 0.00478 0.00478  0.00478
 6 0.01  1     1        0.00478 0.00478 0.00478 0.00478 0.00478 0.00478  0.00478
 7 0.012 1     1        0.00478 0.00478 0.00478 0.00478 0.00478 0.00478  0.00478
 8 0.014 1     1        0.00478 0.00478 0.00478 0.00478 0.00478 0.00478  0.00478
 9 0.016 1     1        0.00478 0.00478 0.00478 0.00478 0.00478 0.00478  0.00478
10 0.018 1     1        0.00478 0.00478 0.00478 0.00478 0.00478 0.00478  0.00478
# ℹ 15,562,490 more rows
# ℹ 1 more variable: `ROC-A2` <dbl>

NOTE: Data frame has 15562500 rows and 11 columns. 

NOTE: Log Print Time:  2023-10-14 13:59:41 
NOTE: Elapsed Time: 15.2787973880768 secs 

Starting artifact detection parameters with step_size = 30 * 5. 

NOTE: Log Print Time:  2023-10-14 13:59:41 
NOTE: Elapsed Time: 0.000774621963500977 secs 

Artifact detection finished. 

NOTE: Log Print Time:  2023-10-14 14:08:13 
NOTE: Elapsed Time: 8.54647764762243 mins 
```
### Automated spindle detection

The library contains implementations of two spindle
detection algorithms discussed in [O'Reilly and Nielsen
(2015)](https://doi.org/10.3389/fnhum.2015.00353). We give a brief overview of
them here but refer to their original publications for further detail.

**Sigma Index Algorithm** : The Sigma Index algorithm [(Huupponen et al.,
2007)](https://pubmed.ncbi.nlm.nih.gov/17555950/) uses the amplitude spectrum to
find spindles by characterizing abnormal values among the spindle frequency
band. Per each $1$ second window of the EEG, $a.$ the
maximum amplitude in the spindle frequency, which we call $S_{max}$, $b.$ the
average amplitude in the low alpha and theta frequencies, which we call
$\alpha_{mean}, \theta_{mean}$, and $c.$ the maximum alpha amplitude
$\alpha_{max}$, are computed. The sigma index is defind to be 

$$f(S_{max}, \alpha_{mean}, \beta_{mean}) = \begin{cases} 
0 & \alpha_{max} > S_{max} \\ 
\frac{2S_{max}}{\alpha_{mean} + \beta_{mean} } & otherwise
\end{cases}$$

Higher values are indicative of a higher spindle probability. The rejection
threshold recommended in the original paper is $\lambda = 4.5$.

To perform spindle detection with the *Sigma index* algorithm, simply call
`eeg$spindle_detection(channel=0, method="sigma_index", filter=TRUE)`

Letting `channel=0` sets spindle detection to be carried out across all
channels. If you want perform it only over the $n$th channel simply let `channel
= n`. If `filter` is `TRUE` then results that do not meet the recommended 
threshold $f \geq 4.5$ are automatically removed.

### Power spectrum analysis

It is straightforward to estimate the power spectral density of the EEG signals
using the package. For the sake of showcasing, we will only show the spectrum of
the first $40$ epochs ($1200$ minutes) of the record.

```r
eeg <- EEG$new('some/relative/path/eeg.csv')
eeg$subset(1, 40) # Subset to epoch interval [1, 40].
eeg$compute_psd() # Computes PSD and saves results in $psd attribute.
eeg$iplot_psd()
```

![](https://i.ibb.co/F5pD0Hn/PSD.png)

(We have artifcially zoomed into this interactive plot so as to display only
frequencies up to 40 Hz.)


## Example artifact detection

We will perform artifact detection and rejection over the first $10$ minutes of
our record. Afterwards, we will plot the spectrograms of both the raw and the
artifact rejected EEGs. Notice that, explanatory comments aside, artifact
rejection is conducted in only two lines.

```r
eeg <- EEG$new('some/relative/path/eeg.csv')
eeg$subset(1, 20)

# 1. Perform artifact detection. Notice that after artifact detection we are 
# filtering the analysis results so as to preserve only those anomalies with 
# a normalized strength greather than 0.4.

eeg$artf_stepwise(step_size = 120)
eeg$sfilter(0.4)

#2. Perform artifact rejection. This creates a new EEG object that is a copy of the 
# original one, except subepochs containing artifacts are removed. 
clean_eeg <- eeg$artf_reject()

# Plot the spectrograms if you wish to compare

clean_eeg$spectrogram(1) # First channel of the artifact rejected EEG. 
eeg$spectrogram(1) # First channel of the original EEG.
```

![](https://i.ibb.co/55PZP0R/comparison.png)

The clean record has a shorter duration, as it is to be expected from the fact
that artifact contaminated epochs were dropped.

