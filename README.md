
# EEG Computational Toolkit
  

> :last_quarter_moon_with_face: Developed at the [Laboratory for the Study of Sleep Slow-wave activity](https://www.med.upenn.edu/slowwavelab/)

A scientific package for computational EEG analysis.

- [EEG Computational Toolkit](#eeg-computational-toolkit)
    - [Installation and import](#installation-and-import)
    - [Loading EEG data](#loading-eeg-data)
    - [EEG Visualization](#eeg-visualization)
    - [Handling EEG data](#handling-eeg-data)
    - [Filtering](#filtering)
    - [Artifact detection](#artifact-detection)
  - [Power spectrum analysis](#power-spectrum-analysis)


### Installation and import


```r
devtools::install_github("slopezpereyra/EEG-toolkit")
library(eegtk)
```
### Loading EEG data

To load EEG data, we use the `load_eeg(data_file, ...)` function. This function takes a data file (with an optional signals file) as arguments.

For example, after transforming a `test.edf` EEG record to `.csv` format, we run

```r
eeg <- load_eeg("test_data.txt",
                "test_signals.txt") %>%
				na.omit()
```

`load_eeg` returns an `eeg` S4 object containing the `@data` and `@signals` atrributes, the latter being an empty data frame if no signals file was provided.
  
### EEG Visualization

We can produce static or interactive plots to get a visual sense of our record. For example, let's take a peek at the thirstiest epoch in our newly created `eeg` object. 

```
epoch <- get_epoch(eeg, 30)
plot(epoch)
```
![enter image description here](https://i.ibb.co/0X4GG8T/Screenshot-from-2022-12-05-13-12-17.png)

We could have also produced an interactive plot to inspect our EEG record, or a portion of it, live. For example,

`iplot(epoch)`


![Alt Text](https://i.ibb.co/0XgxhKv/ezgif-1-2c5fd0d1e6.gif)

  

### Handling EEG data

<<<<<<< HEAD
EEG data is extremely large. For a sampling rate $f_s = 500$, a single $30$ seconds epoch contains $15.000$ observations!  Hence, it is often desirable to work either with subsets of the record, or lower resolution version of it. These package provides functions to make subsetting and resolution-scaling practical and easy.
=======
EEG data is extremely large. For a sampling rate $f_s = 500$, a single $30$ seconds epoch contains $15.000$ observations! An eight-hour record $-$say, from a sleep session$-$ would have $14.400.000$ measures. 

Hence, it is often desirable to work either with subsets of the record, or lower resolution version of it. These package provides functions to make subsetting and resolution-scaling practical and easy.
>>>>>>> main

The `get_epoch` function has already been showcased. But it is nothing more than a wrapper for a specific call of the `subset_eeg(eeg, start, end)` function, where the `start, end` artifacts are numbers referencing time.

Resampling is also easy via `resample(eeg, n)`. Resampling will keep only one every $n$ observations, producing a lower resolution version of the record. This can help accelerate different types of analysis, such as artifact detection and power spectrum analysis, as well as static and interactive plotting.
 
### Filtering
Low-pass, high-pass and bandpass filters are available via the `low_pass, high_pass` and `bandpass` functions. For example,
```r
minute <- subset_eeg(eeg, 60, 120) # Second minute of record
fminute <- low_pass(fminute, 20) # Apply 20Hz filter
```
 
defines two `eeg` objects containing the second minute of record, the latter with a low-pass filter. Visually, they differ as shown below.

![enter image description here](https://i.ibb.co/HnG5jTc/plot-3.png)


### Artifact detection

Artifact detection is carried out via the CAPA statistical method ([Fisch, Eckley & Fearnhead, 2021](https://onlinelibrary.wiley.com/doi/full/10.1002/sam.11586)). CAPA is adapted to the specificities of EEG data via _sub rosa_ operations. To perform artifact detection simply call `artf(eeg, start, end, ...)`. :microscope:
  

The `artf` function takes an `eeg` object and returns an `analysis` object. An `analysis` object has

- a `@canoms` attribute: data frame of all collective anomalies;

- a `@panoms` attribute: data frame of all point anomalies;

- an `@eeg` attribute: containing the eeg upon which the analysis was conducted.

  

An `analysis` object can be plotted to obtain a visual representation of the anomalies found. For example, the thirstiest epoch contains a few unusual spikes:

```
an <- artf(epoch) # Remember epoch <- get_epoch(eeg, 30)
plot(an)
```

![enter image description here](https://i.ibb.co/BjL6fDR/Screenshot-from-2022-12-05-13-18-15.png)

Marked in red are the primary suspects.  The fourth channel is omitted, since no anomaly was found. 

It should be noted that anomalies have different strengths. We can easily filter an analysis object so as to keep only anomalies stronger than a certain threshold value. For example,

```
filtered <- sfilter(an, 0.1) 
# After min-max normalization, 
# keep only anomalies with strength >= 0.1
```

![enter image description here](https://i.ibb.co/djz0v74/Screenshot-from-2022-12-05-13-21-10.png)
  

## Power spectrum analysis

It is straightforward to estimate the power spectral density of the EEG signals using the package. For the sake of showcasing, we will only show the spectrum of the first 1200 minutes of record.

```
s <- subset_eeg(eeg, 0, 1200)
sd <- psd(s) # Estimate the spectral density
iplot_psd(sd).
```

![](https://i.ibb.co/F5pD0Hn/PSD.png)

(We have artifcially zoomed into this interactive plot so as to display only frequencies up to 40 Hz.)

We may also compute the spectogram of a specific EEG channel. For example, here's the spectogram of the C4-A1 (artifact-contaminated) channel from the start of the record to the sixth minute.

```
x <- subset_eeg(eeg, 0, 360)
spectogram(s, channel=4, hcolors=10) 
# hcolors determines number of colors in the palette
```

![enter image description here](https://i.ibb.co/HqJCTDg/Screenshot-from-2022-12-05-13-34-03.png)