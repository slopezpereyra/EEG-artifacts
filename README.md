# EEG-artifacts

*Automated EEG Artifact Detection*

An implementation of Multivariate CAPA to EEG data for automated artifact detection.

### Artifact detection
Though EEG is designed to record cerebral activity, it also records electrical activity arising from sites other than the brain. Non-cerebral EEG measures are termed artifacts. Both physiologic and extra-physiologic artifacts are a form of data contamination that obscure the analysis of brain activity. Though artifact detection is usually done manually by a trained expert, the field is slowly transitioning to automated artifact detection methods.

### CAPA


[CAPA](https://arxiv.org/abs/1806.01947) is a statistical method for the detection of collective and point anomalies. It works by inferring the location of joint epidemic changes in mean and/or variance through the minimization of a suitable cost function. CAPA is a relatively new statistical method and has never before been used on EEG data.