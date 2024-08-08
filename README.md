# ephys_analysis_MATLAB
Set of MATLAB routines for LFP analysis. Sleep architecture, spectral analysis, coherence, cross-frequency coupling, sharp wave-ripple, delta, and spindle detection

# First steps
It is important to have your data already blocked in 10 epochs and the sleep-wake cycle states sorted (check the repository https://github.com/ikaro-beraldo/sleep-stages-classification)
- Create a folder 'blocked_data' and inside it create a folder for each set of data (each day of recording for each subject eg. 'B1_D1', 'B1_D2', 'B2_D1', 'B2_D2') and add your blocked data (as 'blocked_data.mat').
- Add the 'GMM_classification.mat' (sleep-wake cycle classification results) files inside the subjects/days respective folders.

# The first set of functions 
The 'main_routine' function is the main routine.
 - Change the 'folder_names', 'subject_names' cell variables accordingly to your data set.
 - Run it!
 - Obs: it is possible to run all the set of analysis for each subject/day separatelly by defining the _names_ variables as cells with only 1 element (eg. 'B1_D1')

# Spectral analysis 
 - Power Spectral Density
![Figura 1 - Exploring the Signal_Prancheta 1-01](https://github.com/user-attachments/assets/5d045354-89b8-4de9-acd4-ad2a8ddb8add)

# Signal coherence and cross-frequency modulation analysis
 - Magnitude Squared Coherence
 - Phase Coherence
![Figura 1 - Exploring the Signal_Prancheta 1-02](https://github.com/user-attachments/assets/7a8561ce-6cf9-4bd3-82a0-392e9ee7b38d)

# Electrophysiological events
- Sharp-wave ripples
- Delta
- Spindle
![Figura 1 - Exploring the Signal-03](https://github.com/user-attachments/assets/cbf95443-a6f6-4d1f-8cc5-d9ebf9fd4de5)



