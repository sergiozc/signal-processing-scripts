# Overview
This repository contains some simulation macros related to signal processing.
The macros have been built using Matlab and Python.

# Voice Processing
Pitch frequency measurement, periodogram of each vowel, voice spectrum analysis and cepstrum FFT.
This macro can be found at "voice-processing/" directory.

# Kalman filter
Tracking of a geostationary orbit using extended Kalman filtering. A filter is implemented that performs the prediction from a dynamic model and provided measurements.

# Equalizer
The equalization of a QAM-16 signal transmitted through a multipath channel can be found at "equalizer/" directory.
The simulation has been carried out for a static channel (TSE-LS) and for a changing channel (TSE decision-directed).
This simulation can be found at "/equalizer" directory

# Image Processing
At "image-processing/" directory you can find two different macros:
- An implementation of a JPEG compression
- Image enhancement and 2D filters

# DVB-T
This macro implements the decoding of a DVB-T channel, visualizing TPS carriers.
PD: The signal is not provided but it can be recorded by a RTL-SDR
The code can be found at "DVB-T/" directory.

# Beamforming
An acoustic delay and sum beamformer has been implemented in order to enhance an audio signal.
Directivity patterns are also shown.
