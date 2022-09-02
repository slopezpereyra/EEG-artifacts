# Overview

This software is an implementation of [Multivariate-CAPA](https://arxiv.org/abs/1806.01947) for artifact detection on EEG data.

Though EEG is designed to record cerebral activity, it also records electrical activity arising from sites other than the brain. Non-cerebral EEG measures are termed artifacts. Both physiologic and extra-physiologic artifacts are a form of data contamination that obscure the analysis of brain activity. For example, it is deemed scientifically desirable to run power spectral analysis on artifact-free EEG records.

Artifact detection is usually carried out manually by a trained expert, via careful examination of the EEG records. However, the process is slow, requires time and resources, and shows suboptimal inter-rater reliability. As a consequence, automated artifact detection methods are slowly emerging.

Most of this methods rely on ICA and machine learning. However, artifacts can be understood as statistically anomalous segments of data and therefore may be subject to strictly statistical detection —_id est_, detection methods that do not involve learning.

### CAPA

[CAPA](https://arxiv.org/abs/1806.01947) is a statistical method for the detection of collective and point anomalies. It works by inferring the location of joint epidemic changes in mean and/or variance through the minimization of a suitable cost function. CAPA is a fairly new statistical method and has never before been used on EEG data.
