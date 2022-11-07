# EEG Automated Artifact Detection
- [EEG Automated Artifact Detection](#eeg-automated-artifact-detection)
- [1. Requirements](#1-requirements)
- [2. Installation](#2-installation)
- [3. Features](#3-features)
> :last_quarter_moon_with_face: Developed at the [Laboratory for the Study of Sleep Slow-wave activity](https://www.med.upenn.edu/slowwavelab/)

# 1. Requirements

Aside from R, this library depends on Python >= 3.10.8 and the Python libraries `pandas`, `numpy`, `plotly` and `plotly-resampler`. This is due to the fact that interactive plotting is, for techinical reasons, conducted via Python rather than directly on R.

# 2. Installation 

```r
devtools::install_github("slopezpereyra/EEG-artifacts")
```
  
---
# 3. Features

The package provides the EEG tools typically used in neuroscientific research. This includes resampling, filtering, plotting, interactive visualization, artifact detection and power spectrum analysis.

 In particular, artifact detection is carried out via M-CAPA analysis ([Fisch, Eckley & Fearnhead, 2021](https://onlinelibrary.wiley.com/doi/full/10.1002/sam.11586)). Because the package was originally intended only for artifact detection, this is the most covered feature, with a rich set of analysis parameters and several statitc and interactive plotting tools.

  