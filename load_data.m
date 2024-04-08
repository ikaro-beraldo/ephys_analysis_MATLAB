function data = load_data(file)
% Function that load mat data and put it in struct form 
% The data is already windowing in 630 epochs of 19.5s
% 
% Natalia Espinosa 28/06/2023
% ------------------------------------------------------------------------

    data_epochs = load(file);
    data.Fs = 10000; % Sample frequency
    T = 1/data.Fs;   % Sample Period
    size_data  = length(data_epochs.data3(:,1));
    time_epochs = (0:1:size_data-1)*T; % Time stams of data
    
    % Data recorded at Hippocampo
    data.HPC.data_w = data_epochs.data4'; 
    data.HPC.time_w = time_epochs;
    
    % Data recorded at mediam Prefrontal Cortex
    data.mPC.data_w = data_epochs.data3'; 
    data.mPC.time_w = time_epochs;
    
    % Baseline HFS1 and HFS2 mPC
    data.mPC.BL.data_w   = data_epochs.data3(:,1:92)'; 
    data.mPC.BL.time_w   = time_epochs;
    data.mPC.HFS1.data_w = data_epochs.data3(:,93:184)';
    data.mPC.HFS1.time_w = time_epochs;
    data.mPC.HFS2.data_w = data_epochs.data3(:,185:630)';
    data.mPC.HFS2.time_w = time_epochs;

    % Event times
    data.event_time = (1:1:630)*19.5/60; % time in minutes
    HFS1 = find(abs(data.event_time - 30)< 0.3250);
    data.HFS1_time = data.event_time(HFS1(1)); % trial 92
    
    HFS2 = find(abs(data.event_time - 60)< 0.3250);
    data.HFS2_time = data.event_time(HFS2(1)); % trial 184



end