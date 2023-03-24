
# EEG Computational Toolkit
  

> :last_quarter_moon_with_face: Developed at the [Laboratory for the Study of
> Sleep Slow-wave activity](https://www.med.upenn.edu/slowwavelab/)

A scientific package for computational EEG analysis.

- [EEG Computational Toolkit](#eeg-computational-toolkit)
    - [Installation and import](#installation-and-import)
    - [Loading EEG data](#loading-eeg-data)
    - [EEG Visualization](#eeg-visualization)
    - [Resampling and resolution scaling](#resampling-and-resolution-scaling)
    - [Filtering](#filtering)
    - [Artifact detection](#artifact-detection)
  - [Power spectrum analysis](#power-spectrum-analysis)
  - [Example artifact detection](#example-artifact-detection)


### Installation and import


```r 
devtools::install_github("slopezpereyra/EEG-toolkit") 
library(eegtk) 
```

If R complains of not having an up-do-date version of Rtools, you may try to
use the forcing install option.

```r 
devtools::install_github("slopezpereyra/EEG-toolkit", force=TRUE) 
library(eegtk) 
```

### Loading EEG data

To load EEG data, we use the `EEG$new(data_file, ...)` initiation function. 
This function takes a data file (with an optional signals file) as arguments.

For example, if we have an `eeg.csv` file in some path relative to our working
directory, we call

```r 
eeg <- EEG$new("/relative/path/eeg.csv")
```

The function returns an `eeg` R6 object containing the following attributes: 

- `$data` : The EEG data.
- `$signals` : The channels of the EEG, if a signals file was provided on
  initialization. 

- `$canoms` : Collective anomalies (artifacts) found in the EEG. If artifact
  detection was not yet performed, it defaults to an empty data frame.
- `$panoms` : Point anomalies (artifacts) found in the EEG. If artifact
  detection was not yet performed, it defaults to an empty data frame.
- `$psd` : The estimated power spectrum density of the EEG. If PSD computation
  was not yet performed, it defaults to an empty data frame.
- `$fs` : The sampling frequency $f_s$ of the EEG.

Most methods of the EEG class perform inplace replacement. For example, take the
following lines. 

``` r
eeg <- EEG$new("/relative/path/eeg.csv") 
eeg$subset(1, 10)  
```

This subsets the `eeg` object from epoch $1$ to $10$ *inplace*, meaning that the
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
  

## Power spectrum analysis

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
