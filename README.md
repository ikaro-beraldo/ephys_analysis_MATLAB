# ephys_analysis_MATLAB
Set of MATLAB routines for LFP analysis. Sleep architecture, spectral analysis, coherence, cross-frequency coupling, sharp wave-ripple, delta, and spindle detection

The first set of functions (as when Ysgramor settled the first free-man hold in Tamriel.. and taught a lesson to those filthy elves)
- The 'main_routine' is the main routine (lol). You gotta have your data already blocked in 10 epochs and the sleep-wake cycle states sorted. 
- Create a folder 'blocked_data' and inside it create a folder for each set of data (each day of recording for each subject) and add your blocked data (as 'blocked_data.mat')
- Add the 'GMM_classification.mat' files inside the subjects/days respective folders.

Now you change the folder names, subject names, etc.. on the 'main_routine.m' and than run it!
