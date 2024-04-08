%% Delta detection based on Naty's algorithm

% LFP (LFP data divided in blocks (n*m;  n = number of epochs, m = epochs
% length)
% sd_thrshold = peak threshold (default = 2.8)
% plot_fig = define if the detected deltas are going to be plotted
% (default = true)

function [delta_blocks, parameters] = detect_Delta_Naty_algorithm_rem(LFP,fs,sd_threshold,plot_fig)
% Set some parameters
if isempty(sd_threshold)
    sd_threshold = 2.8;     % Peak threshold (peak-onset or peak-offset)
end
if isempty(plot_fig)
    plot_fig = true;        % Define if it is necessary to plot the figures
end
%% LFP

data_windowing.time_w = 1:size(LFP,2);

n = 4; % filt order
nyquist_rate = fs/2;
Wn = [6 8]/nyquist_rate;  % Lower and Upper frequency limits
ftype = 'stop'; % filter type
[b,a] = butter(n,Wn,ftype); % Create the butter filter

% Filtering LFP data
A8_win = filtfilt(b,a,LFP')';

% Filtered data
A8_win = eegfilt2(A8_win,fs,1,4);
% A8_win = eegfilt2(LFP,fs,6,8,0,0,1);

% Z-scored data
data_zscore  = zscore(A8_win,0,2);
% derivated data
data_diff = diff(data_zscore,[],2)/0.001;
% Delta detection counter
k = 1;

% Extract delta frequency, z-score and derivate of data
for i = 1:size(LFP, 1) % i -> Number of trials
    wname = 'sym8'; % wavelet type
    %     [LoD,HiD,LoR,HiR] = wfilters(wname); % Low and High pass filters Decomposition and Reconstruction
    %     [C,L]             = wavedec(LFP(i,:), 8, LoD, HiD); % wavelet decomposition 8 levels
    %     A8_win(i,:)       = wrcoef('a', C, L, wname, 8);  % DELTA
    %data_zscore(i,:)  = zscore(A8_win(i,:));          % z-scored data
    %data_diff(i,:)    = diff(data_zscore(i,:))/0.001; % derivated data

    % Find all zero crossing indeces(indx) on derivate of z-score data
    indx_zerocross = find(data_diff(i,2:end).* data_diff(i,1:end-1)<0 )+1;

    if (isempty(indx_zerocross)==0)
        % Find peak of delta wave
        [~,indx_peak] = findpeaks(data_zscore(i,:));
        [~,~,int_zcross] = intersect(indx_peak,indx_zerocross);

        % ===== Delta detection =====

        % Check delta events for each peak
        for jj = 1:length(int_zcross)

            if ((int_zcross(jj) > 1) && (int_zcross(jj) < length(indx_zerocross)))
                indx_peak_delta   = indx_zerocross(int_zcross(jj));     % Index vector of peak points
                indx_onset_delta  = indx_zerocross(int_zcross(jj) - 1); % Index vector of onset points
                indx_offset_delta = indx_zerocross(int_zcross(jj) + 1); % Index vector of offset points

                % ==== PARAMETERS =======
                Amp_min     = 250; % uV
                Amp_min_sd = sd_threshold; % Standard Deviation
                Amplitude_onset   = abs(A8_win(i,indx_peak_delta) - A8_win(i,indx_onset_delta)); %[uV]
                Amplitude_offset   = abs(A8_win(i,indx_peak_delta) - A8_win(i,indx_offset_delta)); %[uV]
                Amplitude_sd_onset = abs(data_zscore(i,indx_peak_delta) - data_zscore(i,indx_onset_delta)); %[STD]
                Amplitude_sd_offset = abs(data_zscore(i,indx_peak_delta) - data_zscore(i,indx_offset_delta)); %[STD]
                time_peak   = data_windowing.time_w(indx_peak_delta);   % [seconds]
                time_onset  = data_windowing.time_w(indx_onset_delta);  % [seconds]
                time_offset = data_windowing.time_w(indx_offset_delta); % [seconds]
                Duration    = time_offset - time_onset;  % [seconds]

                % Slope onset
                Px = [time_onset; time_peak]; % Poinst of polygon in x
                Py = [A8_win(i,indx_onset_delta); A8_win(i,indx_peak_delta)]; % Poinst of polygon in y
                [coef_asscendet] = polyfit(Px,Py,1);
                Slope_onset      = coef_asscendet(1);

                % Slope offset
                Px_dec = [time_peak; time_offset]; % Poinst of polygon in x
                Py_dec = [A8_win(i,indx_peak_delta); A8_win(i,indx_offset_delta) ]; % Poinst of polygon in y
                [coef_desendent] = polyfit(Px_dec,Py_dec,1);
                Slope_offset     = coef_desendent(1);

                % Area under curve
                x_lim = data_windowing.time_w(indx_onset_delta:indx_offset_delta);
                Area = trapz(x_lim, A8_win(i,indx_onset_delta:indx_offset_delta)); % Area under curve [uV^2]
                Energy = trapz(x_lim, (A8_win(i,indx_onset_delta:indx_offset_delta)).^2); % Signal Energy

                % Well defined delta wave Amplitude minimum of 100 uV and time
                % onset minimum 50ms before event
                %                 if ((Amplitude > Amp_min) && (data_windowing.time_w(indx_onset_delta) > -0.05))
                if (Amplitude_sd_onset > Amp_min_sd) && Amplitude_sd_offset > Amp_min_sd
                    delta_blocks(k,:) = [i, time_onset, time_peak, time_offset];
                    parameters(k,:) = [i, Amplitude_onset, Duration, Slope_onset, abs(Slope_offset), Area, Energy];
                    k = k + 1;
                end

            end
        end
    end
end

%% Plot Delta detected delta waves
% If it is necessary to plot
if plot_fig
    figure
    for i = 1:size(LFP,1)
        subplot(4,1,1)
        plot(LFP(i,:))
        hold on

        % Block delta idx
        delta_idx = find(delta_blocks(:,1) == i);
        for j = 1:length(delta_idx)
            t = delta_blocks(delta_idx(j),2):delta_blocks(delta_idx(j),4);
            plot(t,LFP(delta_blocks(delta_idx(j),1),t))
        end
        ylim([-1000 1000])
        hold off

        subplot(4,1,2)
        plot(A8_win(i,:))
        hold on

        % Block delta idx
        delta_idx = find(delta_blocks(:,1) == i);
        for j = 1:length(delta_idx)
            t = delta_blocks(delta_idx(j),2):delta_blocks(delta_idx(j),4);
            plot(t,A8_win(delta_blocks(delta_idx(j),1),t))
        end
        ylim([-1000 1000])
        hold off

        subplot(4,1,3)
        plot(data_zscore(i,:))
        hold on

        % Block delta idx
        delta_idx = find(delta_blocks(:,1) == i);
        for j = 1:length(delta_idx)
            t = delta_blocks(delta_idx(j),2):delta_blocks(delta_idx(j),4);
            plot(t,data_zscore(delta_blocks(delta_idx(j),1),t))
        end
        ylim([-5 5])
        hold off

        subplot(4,1,4)
        plot(data_diff(i,:))
        hold on

        % Block delta idx
        delta_idx = find(delta_blocks(:,1) == i);
        for j = 1:length(delta_idx)
            t = delta_blocks(delta_idx(j),2):delta_blocks(delta_idx(j),4);
            plot(t,data_diff(delta_blocks(delta_idx(j),1),t))
        end
        ylim([-100 100])
        hold off
        pause

    end
end
end