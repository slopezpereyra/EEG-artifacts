
# EEG Computational Toolkit
  

> :last_quarter_moon_with_face: Developed at the [Laboratory for the Study of Sleep Slow-wave activity](https://www.med.upenn.edu/slowwavelab/)

A scientific package for computational EEG analysis.

- [EEG Computational Toolkit](#eeg-computational-toolkit)
    - [Requirements](#requirements)
    - [Installation](#installation)
    - [Loading EEG data](#loading-eeg-data)
    - [EEG Visualization](#eeg-visualization)
    - [Handling EEG data](#handling-eeg-data)
    - [Artifact detection](#artifact-detection)
  - [Power spectrum analysis](#power-spectrum-analysis)


 ### Requirements

Aside from R, this library depends on Python >= 3.10.8 and the Python libraries `pandas`, `numpy`, `plotly` and `plotly-resampler`. This is due to the fact that interactive plotting is, for techinical reasons, conducted via Python rather than directly on R.

### Installation

  

In order to perform automated artifact detection, we need to install the `artifactor` package in R. Run the following commands in the R console.

  

```r
devtools::install_github("slopezpereyra/EEG-toolkit")
```

The package can be loaded with `library(eegtk)`.
### Loading EEG data

To load EEG data, we use the `load_eeg(data_file, signals_file = NULL)` function. This function takes a data file and an optional signals file as arguments.

For example, after previously exporting from EDF Browser an EDF file named `test`, we run

```r
eeg <- load_eeg("test_data.txt",
                "test_signals.txt") %>%
				na.omit()
```


The `load_eeg` function returns an `eeg` object containing the `@data` and `@signals` atrributes, the latter being an empty data frame if no signals file was provided.

We can consult the sampling frequency if we ignore it.

```
get_sampling_frequency(eeg) # In Hz
> 500 
```

  
### EEG Visualization

We can produce static or interactive plots to get a look at our EEG record. For example, let's take a pick at the thirstiest epoch in our newly created `eeg` object. 

```
epoch <- get_epoch(eeg, 30)
plot(epoch)
```
![enter image description here](https://i.ibb.co/0X4GG8T/Screenshot-from-2022-12-05-13-12-17.png)

Interactive plotting, on the other hand, is performed through the Python distribution of Plotly. Here's a little teaser!


![Alt Text](https://i.ibb.co/0XgxhKv/ezgif-1-2c5fd0d1e6.gif)

  

### Handling EEG data

Due to the size of EEG data, it may be desirable to work with only a subset of the record, or a low resolution version of it. There are many functions useful in reducing the size of the resolution of our data. The `get_epoch` function has already been showcased. But it is nothing more than a wrapper for a specific call of the `subset_eeg(eeg, start, end)` function. 

A resampling function, `resample(eeg, n)`, is also of great use in reducing the resolution of the EEG. Working with lower-resolution versions of our record can help us accelerate different types of analysis, such as artifact detection, as well as static and interactive plotting.
 
We may also apply low-pass, high-pass and bandpass filters. For example,
```r
minute <- subset_eeg(eeg, 60, 120) # Second minute of record
fminute <- low_pass(fminute, 20) # Apply 20Hz filter
```
 
defines two `eeg` objects containing the second minute of record, the latter with a low-pass filter. Visually, they differ as shown in the image below.

![enter image description here](https://i.ibb.co/HnG5jTc/plot-3.png)


### Artifact detection

Artifact detection is carried out via the CAPA statistical method ([Fisch, Eckley & Fearnhead, 2021](https://onlinelibrary.wiley.com/doi/full/10.1002/sam.11586)). CAPA is adapted to the specificities of EEG data via _sub rosa_ operations. To perform artifact detection simply call the `artf(eeg, start, end, ...)`. :microscope:

  

The `artf` function takes an `eeg` object and returns an `analysis` object. An `analysis` object has

  

- a `@canoms` attribute: data frame of all collective anomalies;

- a `@panoms` attribute: data frame of all point anomalies;

- an `@eeg` attribute: containing the eeg upon which the analysis was conducted.

  

The `analysis` class has numerous plotting methods, the most important of which is `plot(analysis)` (used to create the graphs of the forthcoming section).

For example, the thirstiest epoch contains a few unusual spikes:

```
an <- artf(epoch) # Remember epoch <- get_epoch(eeg, 30)
plot(an)
```

![enter image description here](https://i.ibb.co/BjL6fDR/Screenshot-from-2022-12-05-13-18-15.png)

Marked in red are the primary suspects.  The fourth channel is omitted, since no anomaly was found. It should be noted that anomalies have different strengths. We can easily filter an analysis object so as to keep only anomalies stronger than a certain threshold value. For example,

```
filtered <- sfilter(an, 0.1) # After min-max normalization, keep only 
							# anomalies with strength >= 0.1
```

![enter image description here](https://i.ibb.co/djz0v74/Screenshot-from-2022-12-05-13-21-10.png)
  
  It should be added that analysis objects (i.e. artifact detected segments of the EEG) also allow for interactive plotting.

## Power spectrum analysis

It is quite simple to estimate the power spectral density of the EEG signals using the package. For the sake of showcasing, we will only show the spectrum of the first 6 minutes of record.

```
s <- subset_eeg(eeg, 0, 360)
sd <- psd(s) # Estimate the spectral density
plot_psd(sd, xlim=30) # Show only frequencies up to 30 Hz.
```

![enter image description here](https://i.ibb.co/X3TpkQC/Screenshot-from-2022-12-05-13-30-54.png)

We may also compute the spectogram of a specific EEG channel. For example, here's the spectogram of the C4-A1 (artifact-contaminated) channel.

```
spectogram(s, channel=4, hcolors=10) # hcolors determines number of 
									# colors in the palette
```

![enter image description here](https://i.ibb.co/HqJCTDg/Screenshot-from-2022-12-05-13-34-03.png)