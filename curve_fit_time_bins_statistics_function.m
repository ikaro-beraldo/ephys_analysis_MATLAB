%% Statistics function
function [results, normalized_results, anova_results, anova_rejected, filtered_correlation, results_names] = curve_fit_time_bins_statistics_function(directory,names,channel_selection,time_limit,b_length,bin_length)

counter_loop = 0;
results = struct;
if isempty(time_limit)
    time_limit = Inf;
else
    % Time limit to consider in minutes
    time_limit = time_limit * 60 / b_length;
end

% Get each recording length before hand
for ii = 1:size(names,1)
    for jj = 1:size(names,2)

        %% 0 - Get the classification results
        open_class = fullfile(directory,char(names{ii,jj}),'GMM_Classification.mat');
        load(open_class,'GMM','GMM_NREM_All_Sort','GMM_REM_All_Sort','GMM_WK_All_Sort')

        recording_length(ii,jj) = length(GMM.All_Sort);
    end 
end

% Get the maximum length for all recordings
max_rec_length = max(recording_length,[],'all');
clear recording_length

for ii = 1:size(names,1)
    for jj = 1:size(names,2)

        %% 0 - Get the classification results
        open_class = fullfile(directory,char(names{ii,jj}),'GMM_Classification.mat');
        load(open_class,'GMM','GMM_NREM_All_Sort','GMM_REM_All_Sort','GMM_WK_All_Sort')

        % Fix the classification
        NREM = find(GMM.All_Sort == 2);
        REM = find(GMM.All_Sort == 1);
        WK = find(GMM.All_Sort == 3);

        % Filter by time
        NREM = NREM(NREM <= time_limit);
        REM = REM(REM <= time_limit);
        WK = WK(WK <= time_limit);

        %% 1 - Script to extract band power throughout the days

        % Load files
        filename = fullfile(directory,'PSD',char(names{ii,jj}));
        load(filename,'BANDS')

        % mPFC - Delta
        results.bands.LFP1.Delta.NREM(jj,ii,:) = time_bin_average(BANDS.Delta.LFP1, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.bands.LFP1.Delta.REM(jj,ii,:) = time_bin_average(BANDS.Delta.LFP1, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.bands.LFP1.Delta.WK(jj,ii,:) = time_bin_average(BANDS.Delta.LFP1, bin_length, b_length, max_rec_length, WK, true, [], []);

        % mPFC - Theta
        results.bands.LFP1.Theta.NREM(jj,ii,:) = time_bin_average(BANDS.Theta.LFP1, bin_length, b_length, max_rec_length, NREM, true, [], []);        
        results.bands.LFP1.Theta.REM(jj,ii,:) = time_bin_average(BANDS.Theta.LFP1, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.bands.LFP1.Theta.WK(jj,ii,:) = time_bin_average(BANDS.Theta.LFP1, bin_length, b_length, max_rec_length, WK, true, [], []);

        % mPFC - Gamma
        results.bands.LFP1.Gamma.NREM(jj,ii,:) = time_bin_average(BANDS.Gamma.LFP1, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.bands.LFP1.Gamma.REM(jj,ii,:) = time_bin_average(BANDS.Gamma.LFP1, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.bands.LFP1.Gamma.WK(jj,ii,:) = time_bin_average(BANDS.Gamma.LFP1, bin_length, b_length, max_rec_length, WK, true, [], []);

        % CA11 - Delta
        results.bands.LFP2.Delta.NREM(jj,ii,:) = time_bin_average(BANDS.Delta.LFP2, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.bands.LFP2.Delta.REM(jj,ii,:) = time_bin_average(BANDS.Delta.LFP2, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.bands.LFP2.Delta.WK(jj,ii,:) = time_bin_average(BANDS.Delta.LFP2, bin_length, b_length, max_rec_length, WK, true, [], []);

        % CA11 - Theta
        results.bands.LFP2.Theta.NREM(jj,ii,:) = time_bin_average(BANDS.Theta.LFP2, bin_length, b_length, max_rec_length, NREM, true, [], []);        
        results.bands.LFP2.Theta.REM(jj,ii,:) = time_bin_average(BANDS.Theta.LFP2, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.bands.LFP2.Theta.WK(jj,ii,:) = time_bin_average(BANDS.Theta.LFP2, bin_length, b_length, max_rec_length, WK, true, [], []);

        % CA11 - Gamma
        results.bands.LFP2.Gamma.NREM(jj,ii,:) = time_bin_average(BANDS.Gamma.LFP2, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.bands.LFP2.Gamma.REM(jj,ii,:) = time_bin_average(BANDS.Gamma.LFP2, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.bands.LFP2.Gamma.WK(jj,ii,:) = time_bin_average(BANDS.Gamma.LFP2, bin_length, b_length, max_rec_length, WK, true, [], []);

        % CA12 - Delta
        results.bands.LFP3.Delta.NREM(jj,ii,:) = time_bin_average(BANDS.Delta.LFP3, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.bands.LFP3.Delta.REM(jj,ii,:) = time_bin_average(BANDS.Delta.LFP3, bin_length, b_length,  max_rec_length, REM, true, [], []);
        results.bands.LFP3.Delta.WK(jj,ii,:) = time_bin_average(BANDS.Delta.LFP3, bin_length, b_length,  max_rec_length, WK, true, [], []);

        % CA12 - Theta
        results.bands.LFP3.Theta.NREM(jj,ii,:) = time_bin_average(BANDS.Theta.LFP3, bin_length, b_length, max_rec_length, NREM, true, [], []);        
        results.bands.LFP3.Theta.REM(jj,ii,:) = time_bin_average(BANDS.Theta.LFP3, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.bands.LFP3.Theta.WK(jj,ii,:) = time_bin_average(BANDS.Theta.LFP3, bin_length, b_length, max_rec_length, WK, true, [], []);

        % CA12 - Gamma
        results.bands.LFP3.Gamma.NREM(jj,ii,:) = time_bin_average(BANDS.Gamma.LFP3, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.bands.LFP3.Gamma.REM(jj,ii,:) = time_bin_average(BANDS.Gamma.LFP3, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.bands.LFP3.Gamma.WK(jj,ii,:) = time_bin_average(BANDS.Gamma.LFP3, bin_length, b_length, max_rec_length, WK, true, [], []);

        %% 2 - Script to extract phase coherence throughout the days

        % Load files
        if channel_selection(ii) == 2
            filename = fullfile(directory,'phase_coherence12',char(names{ii,jj}));
        elseif channel_selection(ii) == 3
            filename = fullfile(directory,'phase_coherence13',char(names{ii,jj}));
        end

        load(filename,'Bands')

        % Delta
        results.phase.Delta.NREM(jj,ii,:) = time_bin_average(Bands.Delta, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.phase.Delta.REM(jj,ii,:) = time_bin_average(Bands.Delta, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.phase.Delta.WK(jj,ii,:) = time_bin_average(Bands.Delta, bin_length, b_length, max_rec_length, WK, true, [], []);

        % Theta
        results.phase.Theta.NREM(jj,ii,:) = time_bin_average(Bands.Theta, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.phase.Theta.REM(jj,ii,:) = time_bin_average(Bands.Theta, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.phase.Theta.WK(jj,ii,:) = time_bin_average(Bands.Theta, bin_length, b_length, max_rec_length, WK, true, [], []);

        % Gamma
        results.phase.Gamma.NREM(jj,ii,:) = time_bin_average(Bands.Gamma, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.phase.Gamma.REM(jj,ii,:) = time_bin_average(Bands.Gamma, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.phase.Gamma.WK(jj,ii,:) = time_bin_average(Bands.Gamma, bin_length, b_length, max_rec_length, WK, true, [], []);


        %% 3 - Script to extract spectral coherence throughout the days

        % Load files]
        if channel_selection(ii) == 2
            filename = fullfile(directory,'spectral_coherence12',char(names{ii,jj}));
            load(filename,'C_band12')
            C_band = C_band12;
        elseif channel_selection(ii) == 3
            filename = fullfile(directory,'spectral_coherence13',char(names{ii,jj}));
            load(filename,'C_band13')
            C_band = C_band13;
        end

        % Delta
        results.spectral.Delta.NREM(jj,ii,:) = time_bin_average(C_band.Delta, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.spectral.Delta.REM(jj,ii,:) = time_bin_average(C_band.Delta, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.spectral.Delta.WK(jj,ii,:) = time_bin_average(C_band.Delta, bin_length, b_length, max_rec_length, WK, true, [], []);

        % Theta
        results.spectral.Theta.NREM(jj,ii,:) = time_bin_average(C_band.Theta, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.spectral.Theta.REM(jj,ii,:) = time_bin_average(C_band.Theta, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.spectral.Theta.WK(jj,ii,:) = time_bin_average(C_band.Theta, bin_length, b_length, max_rec_length, WK, true, [], []);

        % Gamma
        results.spectral.Gamma.NREM(jj,ii,:) = time_bin_average(C_band.Gamma, bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.spectral.Gamma.REM(jj,ii,:) = time_bin_average(C_band.Gamma, bin_length, b_length, max_rec_length, REM, true, [], []);
        results.spectral.Gamma.WK(jj,ii,:) = time_bin_average(C_band.Gamma, bin_length, b_length, max_rec_length, WK, true, [], []);


        %% 4 - Script to extract MI throughout the days

        % Load files
        filename = fullfile(directory,'MI',char(names{ii,jj}));
        load(filename,'MI_all')

        Pf = [5 12];
        % Amplitude
        Af = [30 55; 65 90; 90 140; 150 170; 180 200; 200 250];

        % LFP2 - 30 - 55
        results.MI.LFP2.f_30_55.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(1,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP2.f_30_55.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(1,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP2.f_30_55.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(1,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP2 - 65 - 90
        results.MI.LFP2.f_65_90.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(2,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP2.f_65_90.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(2,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP2.f_65_90.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(2,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP2 - 90 - 140
        results.MI.LFP2.f_90_140.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(3,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP2.f_90_140.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(3,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP2.f_90_140.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(3,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP2 - 150 - 170
        results.MI.LFP2.f_150_170.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(4,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP2.f_150_170.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(4,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP2.f_150_170.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(4,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP2 - 180 - 200
        results.MI.LFP2.f_180_200.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(5,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP2.f_180_200.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(5,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP2.f_180_200.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(5,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP2 - 200 - 250
        results.MI.LFP2.f_200_250.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(6,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP2.f_200_250.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(6,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP2.f_200_250.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(6,:), bin_length, b_length, max_rec_length, WK, true, [], []);


        % LFP3 - 30 - 55
        results.MI.LFP3.f_30_55.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(1,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP3.f_30_55.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(1,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP3.f_30_55.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(1,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP3 - 65 - 90
        results.MI.LFP3.f_65_90.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(2,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP3.f_65_90.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(2,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP3.f_65_90.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(2,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP3 - 90 - 140
        results.MI.LFP3.f_90_140.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(3,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP3.f_90_140.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(3,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP3.f_90_140.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(3,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP3 - 150 - 170
        results.MI.LFP3.f_150_170.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(4,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP3.f_150_170.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(4,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP3.f_150_170.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(4,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP3 - 180 - 200
        results.MI.LFP3.f_180_200.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(5,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP3.f_180_200.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(5,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP3.f_180_200.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(5,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP3 - 200 - 250
        results.MI.LFP3.f_200_250.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(6,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI.LFP3.f_200_250.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(6,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI.LFP3.f_200_250.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(6,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        %% Delta MI throught days

        % Load files
        filename = fullfile(directory,'MI_delta',char(names{ii,jj}));
        load(filename,'MI_all')

        Pf = [5 12];
        % Amplitude
        Af = [30 55; 65 90; 90 140; 150 170; 180 200; 200 250];

        % LFP2 - 30 - 55
        results.MI_delta.LFP2.f_30_55.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(1,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP2.f_30_55.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(1,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP2.f_30_55.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(1,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP2 - 65 - 90
        results.MI_delta.LFP2.f_65_90.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(2,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP2.f_65_90.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(2,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP2.f_65_90.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(2,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP2 - 90 - 140
        results.MI_delta.LFP2.f_90_140.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(3,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP2.f_90_140.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(3,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP2.f_90_140.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(3,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP2 - 150 - 170
        results.MI_delta.LFP2.f_150_170.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(4,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP2.f_150_170.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(4,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP2.f_150_170.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(4,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP2 - 180 - 200
        results.MI_delta.LFP2.f_180_200.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(5,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP2.f_180_200.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(5,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP2.f_180_200.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(5,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP2 - 200 - 250
        results.MI_delta.LFP2.f_200_250.NREM(jj,ii,:) = time_bin_average(MI_all.LFP2(6,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP2.f_200_250.REM(jj,ii,:) = time_bin_average(MI_all.LFP2(6,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP2.f_200_250.WK(jj,ii,:) = time_bin_average(MI_all.LFP2(6,:), bin_length, b_length, max_rec_length, WK, true, [], []);


        % LFP3 - 30 - 55
        results.MI_delta.LFP3.f_30_55.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(1,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP3.f_30_55.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(1,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP3.f_30_55.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(1,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP3 - 65 - 90
        results.MI_delta.LFP3.f_65_90.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(2,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP3.f_65_90.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(2,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP3.f_65_90.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(2,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP3 - 90 - 140
        results.MI_delta.LFP3.f_90_140.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(3,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP3.f_90_140.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(3,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP3.f_90_140.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(3,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP3 - 150 - 170
        results.MI_delta.LFP3.f_150_170.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(4,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP3.f_150_170.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(4,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP3.f_150_170.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(4,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP3 - 180 - 200
        results.MI_delta.LFP3.f_180_200.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(5,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP3.f_180_200.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(5,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP3.f_180_200.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(5,:), bin_length, b_length, max_rec_length, WK, true, [], []);

        % LFP3 - 200 - 250
        results.MI_delta.LFP3.f_200_250.NREM(jj,ii,:) = time_bin_average(MI_all.LFP3(6,:), bin_length, b_length, max_rec_length, NREM, true, [], []);
        results.MI_delta.LFP3.f_200_250.REM(jj,ii,:) = time_bin_average(MI_all.LFP3(6,:), bin_length, b_length, max_rec_length, REM, true, [], []);
        results.MI_delta.LFP3.f_200_250.WK(jj,ii,:) = time_bin_average(MI_all.LFP3(6,:), bin_length, b_length, max_rec_length, WK, true, [], []);


        %% 5 - Script to extract RDS throughout the days

        % Load files
        filename = fullfile(directory,'ripple',char(names{ii,jj}));
        filename2 = fullfile(directory,'RDS',char(names{ii,jj}));
        if channel_selection(ii) == 2
            load(filename,'ripple_result_block_2')
            load(filename2,'block_2')
            ripple_block = ripple_result_block_2(:,1);
            RDS_block.delta_ripple = block_2.delta_ripple_timestamps.delta(:,1);
            RDS_block.ripple_delta = block_2.ripple_delta_timestamps.ripple(:,1);
            RDS_block.delta_spindle = block_2.delta_spindle_timestamps.delta(:,1);
            RDS_block.ripple_delta_spindle = block_2.ripple_delta_spindle_timestamps.ripple(:,1);
        elseif channel_selection(ii) == 3
            load(filename,'ripple_result_block_3')
            load(filename2,'block_3')
            ripple_blocks = ripple_result_block_3(:,1);
            RDS_blocks.delta_ripple = block_3.delta_ripple_timestamps.delta(:,1);
            RDS_blocks.ripple_delta = block_3.ripple_delta_timestamps.ripple(:,1);
            RDS_blocks.delta_spindle = block_3.delta_spindle_timestamps.delta(:,1);
            RDS_blocks.ripple_delta_spindle = block_3.ripple_delta_spindle_timestamps.ripple(:,1);
        end

        % Delta
        filename = fullfile(directory,'delta',char(names{ii,jj}));
        load(filename,'delta_blocks')
        delta_blocks = delta_blocks(:,1);
        % Spindle
        filename = fullfile(directory,'spindle',char(names{ii,jj}));
        load(filename,'spindle_blocks')
        spindle_blocks = spindle_blocks(:,1);

        % CA1-1 Delta
        filename = fullfile(directory,'delta_CA1',char(names{ii,jj}));
        load(filename,'delta_blocks_2')
        delta_blocks_2 = delta_blocks_2(:,1);

        % CA1-2 Delta
        filename = fullfile(directory,'delta_CA1',char(names{ii,jj}));
        load(filename,'delta_blocks_3')
        delta_blocks_3 = delta_blocks_3(:,1);

        n_epochs = length(GMM.All_Sort);
        results.RDS.ripple(jj,ii,:) = 60*time_bin_counting(ripple_blocks(ripple_blocks <= time_limit), bin_length, b_length, max_rec_length, n_epochs , NREM, true, true);
        results.RDS.delta(jj,ii,:) = 60*time_bin_counting(delta_blocks(delta_blocks <= time_limit), bin_length, b_length, max_rec_length, n_epochs , NREM, true, true);
        results.RDS.spindle(jj,ii,:) = 60*time_bin_counting(spindle_blocks(spindle_blocks <= time_limit), bin_length, b_length, max_rec_length, n_epochs , NREM, true, true);
        results.RDS.ripple_delta(jj,ii,:) = 60*time_bin_counting(RDS_blocks.ripple_delta(RDS_blocks.ripple_delta <= time_limit), bin_length, b_length, max_rec_length, n_epochs , NREM, true, true);
        results.RDS.delta_ripple(jj,ii,:) = 60*time_bin_counting(RDS_blocks.delta_ripple(RDS_blocks.delta_ripple <= time_limit), bin_length, b_length, max_rec_length, n_epochs , NREM, true, true);
        results.RDS.delta_spindle(jj,ii,:) = 60*time_bin_counting(RDS_blocks.delta_spindle(RDS_blocks.delta_spindle <= time_limit), bin_length, b_length, max_rec_length, n_epochs , NREM, true, true);
        results.RDS.ripple_delta_spindle(jj,ii,:) = 60*time_bin_counting(RDS_blocks.ripple_delta_spindle(RDS_blocks.ripple_delta_spindle <= time_limit), bin_length, b_length, max_rec_length, n_epochs , NREM, true, true);
        results.RDS.delta_CA11(jj,ii,:) = 60*time_bin_counting(delta_blocks_2(delta_blocks_2 <= time_limit), bin_length, b_length, max_rec_length, n_epochs , NREM, true, true);
        results.RDS.delta_CA12(jj,ii,:) = 60*time_bin_counting(delta_blocks_3(delta_blocks_3 <= time_limit), bin_length, b_length, max_rec_length, n_epochs , NREM, true, true);

        %% Extract RDS parameters

        % RIPPLE
        filename = fullfile(directory,'ripple',char(names{ii,jj}));
        load(filename,'ripple_result_block_2','ripple_result_block_3','duration','RMS','PEAK','FREQ')
        time_select_2 = find(ripple_result_block_2(:,1) <= time_limit);
        time_select_3 = find(ripple_result_block_3(:,1) <= time_limit);

        results.RDS.duration.ripple_2(jj,ii,:) = time_bin_average(duration.ripple_2(time_select_2), bin_length, b_length, max_rec_length, NREM, true, ripple_result_block_2(:,1), length(GMM.All_Sort));
        results.RDS.duration.ripple_3(jj,ii,:) = time_bin_average(duration.ripple_3(time_select_3), bin_length, b_length, max_rec_length, NREM, true, ripple_result_block_3(:,1), length(GMM.All_Sort));
        results.RDS.RMS.ripple_2(jj,ii,:) = time_bin_average(RMS.ripple_2(time_select_2), bin_length, b_length, max_rec_length, NREM, true, ripple_result_block_2(:,1), length(GMM.All_Sort));
        results.RDS.RMS.ripple_3(jj,ii,:) = time_bin_average(RMS.ripple_3(time_select_3), bin_length, b_length, max_rec_length, NREM, true, ripple_result_block_3(:,1), length(GMM.All_Sort));
        results.RDS.PEAK.ripple_2(jj,ii,:) = time_bin_average(PEAK.ripple_2(time_select_2), bin_length, b_length, max_rec_length, NREM, true, ripple_result_block_2(:,1), length(GMM.All_Sort));
        results.RDS.PEAK.ripple_3(jj,ii,:) = time_bin_average(PEAK.ripple_3(time_select_3), bin_length, b_length, max_rec_length, NREM, true, ripple_result_block_3(:,1), length(GMM.All_Sort));
        results.RDS.FREQ.ripple_2(jj,ii,:) = time_bin_average(FREQ.ripple_2(time_select_2), bin_length, b_length, max_rec_length, NREM, true, ripple_result_block_2(:,1), length(GMM.All_Sort));
        results.RDS.FREQ.ripple_3(jj,ii,:) = time_bin_average(FREQ.ripple_3(time_select_3), bin_length, b_length, max_rec_length, NREM, true, ripple_result_block_3(:,1), length(GMM.All_Sort));

        % SPINDLE
        filename = fullfile(directory,'spindle',char(names{ii,jj}));
        load(filename,'spindle_blocks','duration','RMS','PEAK','FREQ')
        time_select_2 = find(spindle_blocks(:,1) <= time_limit);

        results.RDS.duration.spindle(jj,ii,:) = time_bin_average(duration.spindle(time_select_2), bin_length, b_length, max_rec_length, NREM, true, spindle_blocks(:,1), length(GMM.All_Sort));
        results.RDS.RMS.spindle(jj,ii,:) = time_bin_average(RMS.spindle(time_select_2), bin_length, b_length, max_rec_length, NREM, true, spindle_blocks(:,1), length(GMM.All_Sort));
        results.RDS.PEAK.spindle(jj,ii,:) = time_bin_average(PEAK.spindle(time_select_2), bin_length, b_length, max_rec_length, NREM, true, spindle_blocks(:,1), length(GMM.All_Sort));
        results.RDS.FREQ.spindle(jj,ii,:) = time_bin_average(FREQ.spindle(time_select_2), bin_length, b_length, max_rec_length, NREM, true, spindle_blocks(:,1), length(GMM.All_Sort));

        % DELTA
        filename = fullfile(directory,'delta',char(names{ii,jj}));
        load(filename,'delta_blocks','delta_parameters')
        time_select_2 = find(delta_blocks(:,1) <= time_limit);

        results.RDS.Amplitude_onset.delta(jj,ii,:) = time_bin_average(delta_parameters.parameters(time_select_2,2), bin_length, b_length, max_rec_length, NREM, true, delta_blocks(:,1), length(GMM.All_Sort));
        results.RDS.duration.delta(jj,ii,:) = time_bin_average(delta_parameters.parameters(time_select_2,3), bin_length, b_length, max_rec_length, NREM, true, delta_blocks(:,1), length(GMM.All_Sort));
        results.RDS.Slope_onset.delta(jj,ii,:) = time_bin_average(delta_parameters.parameters(time_select_2,4), bin_length, b_length, max_rec_length, NREM, true, delta_blocks(:,1), length(GMM.All_Sort));
        results.RDS.Slope_offset.delta(jj,ii,:) = time_bin_average(delta_parameters.parameters(time_select_2,5), bin_length, b_length, max_rec_length, NREM, true, delta_blocks(:,1), length(GMM.All_Sort));
        results.RDS.Area.delta(jj,ii,:) = time_bin_average(delta_parameters.parameters(time_select_2,6), bin_length, b_length, max_rec_length, NREM, true, delta_blocks(:,1), length(GMM.All_Sort));
        results.RDS.Energy.delta(jj,ii,:) = time_bin_average(delta_parameters.parameters(time_select_2,7), bin_length, b_length, max_rec_length, NREM, true, delta_blocks(:,1), length(GMM.All_Sort));
        % DELTA CA1
        filename = fullfile(directory,'delta_CA1',char(names{ii,jj}));
        load(filename,'delta_blocks_2','delta_parameters_2','delta_blocks_3','delta_parameters_3')
        time_select_2 = find(delta_blocks_2(:,1) <= time_limit);
        time_select_3 = find(delta_blocks_3(:,1) <= time_limit);

        % Delta CA11
        results.RDS.Amplitude_onset.delta_CA11(jj,ii,:) = time_bin_average(delta_parameters_2.parameters(time_select_2,2), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_2(:,1), length(GMM.All_Sort));
        results.RDS.duration.delta_CA11(jj,ii,:) = time_bin_average(delta_parameters_2.parameters(time_select_2,3), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_2(:,1), length(GMM.All_Sort));
        results.RDS.Slope_onset.delta_CA11(jj,ii,:) = time_bin_average(delta_parameters_2.parameters(time_select_2,4), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_2(:,1), length(GMM.All_Sort));
        results.RDS.Slope_offset.delta_CA11(jj,ii,:) = time_bin_average(delta_parameters_2.parameters(time_select_2,5), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_2(:,1), length(GMM.All_Sort));
        results.RDS.Area.delta_CA11(jj,ii,:) = time_bin_average(delta_parameters_2.parameters(time_select_2,6), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_2(:,1), length(GMM.All_Sort));
        results.RDS.Energy.delta_CA11(jj,ii,:) = time_bin_average(delta_parameters_2.parameters(time_select_2,7), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_2(:,1), length(GMM.All_Sort));

        % Delta CA12
        results.RDS.Amplitude_onset.delta_CA12(jj,ii,:) = time_bin_average(delta_parameters_3.parameters(time_select_3,2), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_3(:,1), length(GMM.All_Sort));
        results.RDS.duration.delta_CA12(jj,ii,:) = time_bin_average(delta_parameters_3.parameters(time_select_3,3), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_3(:,1), length(GMM.All_Sort));
        results.RDS.Slope_onset.delta_CA12(jj,ii,:) = time_bin_average(delta_parameters_3.parameters(time_select_3,4), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_3(:,1), length(GMM.All_Sort));
        results.RDS.Slope_offset.delta_CA12(jj,ii,:) = time_bin_average(delta_parameters_3.parameters(time_select_3,5), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_3(:,1), length(GMM.All_Sort));
        results.RDS.Area.delta_CA12(jj,ii,:) = time_bin_average(delta_parameters_3.parameters(time_select_3,6), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_3(:,1), length(GMM.All_Sort));
        results.RDS.Energy.delta_CA12(jj,ii,:) = time_bin_average(delta_parameters_3.parameters(time_select_3,7), bin_length, b_length, max_rec_length, NREM, true, delta_blocks_3(:,1), length(GMM.All_Sort));

        

        %% 5 - Script to extract RDS throughout the days
        %
        %         % Load files
        %         filename = fullfile(directory,'RDS','totalRDS');
        %         load(filename,'relative')
        %
        %         if channel_selection(ii) == 3
        %             results.RDS.ripple(jj,ii) = relative.ripple_3(ii,jj);
        %             results.RDS.delta(jj,ii) = relative.delta(ii,jj);
        %             results.RDS.spindle(jj,ii) = relative.spindle(ii,jj);
        %             results.RDS.ripple_delta(jj,ii) = relative.ripple_delta_3(ii,jj);
        %             results.RDS.delta_ripple(jj,ii) = relative.delta_ripple_3(ii,jj);
        %             results.RDS.ripple_delta_spindle(jj,ii) = relative.ripple_delta_spindle_3(ii,jj);
        %         elseif channel_selection(ii) == 2
        %             results.RDS.ripple(jj,ii) = relative.ripple_2(ii,jj);
        %             results.RDS.delta(jj,ii) = relative.delta(ii,jj);
        %             results.RDS.spindle(jj,ii) = relative.spindle(ii,jj);
        %             results.RDS.ripple_delta(jj,ii) = relative.ripple_delta_2(ii,jj);
        %             results.RDS.delta_ripple(jj,ii) = relative.delta_ripple_2(ii,jj);
        %             results.RDS.ripple_delta_spindle(jj,ii) = relative.ripple_delta_spindle_2(ii,jj);
        %         end
        %

        %% 6 - Script to extract Sleep Parameters throughout the days

        % Load files
        filename = fullfile(directory,'sleep_parameters',char(names{ii,jj}));
        load(filename,'sleep_param')

        % mPFC - Delta
        results.sleep_parameters.Sleep_E(jj,ii) = sleep_param.Sleep_E;
        results.sleep_parameters.latency_sleep(jj,ii) = sleep_param.latency_sleep1;
        results.sleep_parameters.latency_REM1(jj,ii) = sleep_param.latency_REM1;
        results.sleep_parameters.Sleep_wake_Total(jj,ii) = sleep_param.Sleep_wake_Total;
        results.sleep_parameters.NREM_wake(jj,ii) = sleep_param.NREM_wake;
        results.sleep_parameters.REM_wake(jj,ii) = sleep_param.REM_wake;

        %% 7 - Script to extract architecture throughout the days

        % Load files
        filename = fullfile(directory,char(names{ii,jj}),'architecture');
        load(filename)

        % Delta
        results.architecture.time_spent.NREM(jj,ii) = Time_Spent_NREM;
        results.architecture.time_spent.REM(jj,ii) = Time_Spent_REM;
        results.architecture.time_spent.WK(jj,ii) = Time_Spent_AWAKE;

        % Theta
        results.architecture.n_bouts.NREM(jj,ii) = Number_Bouts_NREM;
        results.architecture.n_bouts.REM(jj,ii) = Number_Bouts_REM;
        results.architecture.n_bouts.WK(jj,ii) = Number_Bouts_AWAKE;

        % Gamma
        results.architecture.d_bouts.NREM(jj,ii) = Duration_Bouts_NREM;
        results.architecture.d_bouts.REM(jj,ii) = Duration_Bouts_REM;
        results.architecture.d_bouts.WK(jj,ii) = Duration_Bouts_AWAKE;

    end
end

%% Get the Barnes Maze results
load("E:\\Barnes Maze - Mestrad\\Resultados DLC - Barnes\\barnes_maze_result_matrix.mat")
% Add to the result struct

results.barnes.distance_trial_av = fix_0_values_on_1st_day(distance_trial_av);
results.barnes.latency_trial_av = fix_0_values_on_1st_day(latency_trial_av);
results.barnes.p_err_trial_av = p_err_trial_sum;
results.barnes.s_err_trial_av = s_err_trial_sum;
results.barnes.speed_trial_av = fix_0_values_on_1st_day(speed_trial_av);
results.barnes.time_on_target_trial_av = fix_0_values_on_1st_day(time_on_target_trial_av);
results.barnes.spatial_trial_av = spatial_trial_sum;
results.barnes.random_trial_av = random_trial_sum;
results.barnes.serial_trial_av = serial_trial_sum;

%% Normalize data by the first day or
% Normalize for the first hour of the first day

% Get every single results field
results_names = fieldnamesr(results,'prefix');

for ii = 1:length(results_names)
    cm = sprintf('matrix = %s;',results_names{ii});
    % execute command
    eval(cm)
    
    if ndims(matrix) == 3
    % Normalize for the first hour of the first day
    matrix = matrix./matrix(1,:,1);
    elseif ndims(matrix) == 2
    % Normalize for the first day
    matrix = matrix./matrix(1,:);
    end
    cm = sprintf('normalized_%s = matrix;', results_names{ii});
    eval(cm)

    % Normalize zscore
    %     matrix = zscore(matrix,[],1);
    %     cm = sprintf('normalized_%s = matrix;', results_names{ii});
    %     eval(cm)
end

%% Calculate anova

% Get every single results field
results_names = fieldnamesr(normalized_results,'prefix');
% Create a table for anova results
anova_results = cell(length(results_names),1);

hour_limit = 6;

% Get a reference group
eval(sprintf('results_reference = %s(:,:,1:hour_limit);',results_names{1}))
% Day group
day = repmat(1:size(names,1),1,numel(results_reference)/size(names,1));
% Animal group
animal = 0;
for ii=1:size(names,2)
    animal = [animal ii*ones(1,size(names,1))];
end
animal = repmat(animal,1,size(results_reference,3));
% Hour group
hour = [];
for ii=1:size(results_reference,3)
    hour = [hour ii*ones(1,numel(names))];
end
% Define the anova groups
groups = {day,hour};
% data = reshape(normalized_results.RDS.ripple(:,:,1:hour_limit),1,[]);
% [p,tbl,stats] = anovan(data,groups,'model','full','varnames',{'Days','Hours'});

for ii = 1:length(results_names)

    % Check if it is going to be a one way or two way anova
    eval(sprintf('check_result = %s;',results_names{ii}))
    % If it is two way
    if ndims(check_result) == 3
        eval(sprintf('data = reshape(%s(:,:,1:hour_limit),1,[]);',results_names{ii}))
        [anova_results{ii,1},anova_results{ii,2},anova_results{ii,3}] = anovan(data,groups,'display','off','model','full','varnames',{'Days','Hours'});
    else % If it is one way
        % execute command
        eval(sprintf("[anova_results{ii,1},anova_results{ii,2},anova_results{ii,3}] = anova1(%s',[],'off');",results_names{ii}))
    end
end

% Get the results names
anova_results = [anova_results results_names];

% Get the comparisons which wielded a p-value < 0.05
rejected = [];
for ii = 1:size(anova_results,1)
    if sum(anova_results{ii,1} < 0.05) > 0  % Check whether one of the anova factors presented a significate difference
        rejected = [rejected; ii];
    end
end
anova_rejected = {results_names{rejected}; anova_results{rejected,1}; anova_results{rejected,2}; anova_results{rejected,3}}';

%% Calculate Pearson Correlation for each combination

% Get every single results field
results_names = fieldnamesr(results,'prefix');
correlation_results = cell(length(results_names)*length(results_names),4);
counter = 1;
for ii = 1:length(results_names)
    for jj = 1:length(results_names)
        eval(sprintf('a = %s;',results_names{ii}));
        eval(sprintf('b = %s;',results_names{jj}));

        % Check whether this is a comparison between a 3d and a 2d arrays
        if ndims(a) ~= ndims(b)
            continue
        end

        % Check if a and b have the same length
        if size(a,1) ~= size(b,1)
            if size(a,1) > size(b,1)
                a = a(1:4,:);
            else
                b = b(1:4,:);
            end
        end

        %         a = nanmean(a,2);
        %         b = nanmean(b,2);
        a = reshape(a,[],1);
        b = reshape(b,[],1);

        [correlation_results{counter,1},correlation_results{counter,2}] = corrcoef(a,b);
        % Fix
        correlation_results{counter,1} = correlation_results{counter,1}(2);
        correlation_results{counter,2} = correlation_results{counter,2}(2);
        correlation_results{counter,3} = results_names{ii};
        correlation_results{counter,4} = results_names{jj};

        counter = counter + 1;
    end
end


%% Bonferroni Correction
n_tests = length([correlation_results{:,2}]);
for ii = 1:n_tests
    correlation_results{ii,2} = correlation_results{ii,2}*n_tests;
end

%% Filter correlation results

clear filtered_correlation
alfa = 0.05;
filtered = find([correlation_results{:,2}] < alfa & [correlation_results{:,2}] ~= 0 & [correlation_results{:,1}] ~= 1);

for ii = 1:length(filtered)
    filtered_correlation{ii,1} = correlation_results{filtered(ii),1};
    filtered_correlation{ii,2} = correlation_results{filtered(ii),2};
    filtered_correlation{ii,3} = correlation_results{filtered(ii),3};
    filtered_correlation{ii,4} = correlation_results{filtered(ii),4};
end
end

function corrected_matrix = fix_0_values_on_1st_day(matrix)
% When the first day has any value == 0, change it to 0.0000001 so it can
% be normalized
matrix(1,matrix(1,:) == 0) = 0.1;
% Corrected matrix
corrected_matrix = matrix;
end

%% Function to average the temporal series in time bins
function average = time_bin_average(data, bin_length, block_length, max_l_expected, state_classification, consider_rest, counts, n_epochs)

% Check whether is necessary to consider the resting elements (so the last
% bin has a lower number of epochs when compared to the rest)
if isempty(consider_rest)
    consider_rest = false;
end

% Whether the total of epochs has not been defined (use the data length)
if isempty(n_epochs)
    l_data = length(data);
else
    l_data = n_epochs;
end

% Define the beginning and end of each bin
b_length_block = (bin_length*3600)/block_length;
% The last bin considers a lower number of epochs
if consider_rest
    n_bins = ceil(l_data/b_length_block);
    beg = ((1:n_bins)-1)*b_length_block + 1;
    en = [(1:n_bins-1)*b_length_block l_data];
    limits = [beg' en'];
% The uneven bin is excluded, so all of them have the same length
else
    n_bins = floor(l_data/b_length_block);
    beg = ((1:n_bins)-1)*b_length_block + 1;
    en = (1:n_bins)*b_length_block;
    limits = [beg' en'];
end


% If the maximum length has been defined
if ~isempty(max_l_expected)
    if consider_rest
        average = nan(1,ceil(max_l_expected/b_length_block));
    else
        average = nan(1,floor(max_l_expected/b_length_block));
    end
else % If there is not an maximum epoch number expected
    average = nan(1,n_bins);
end


% Average stuff
for i = 1:n_bins
    % Firstly check whether is necessary to separate a specific state
    if ~isempty(state_classification)
        bin_values = intersect(limits(i,1):limits(i,2),state_classification);
    else % If it is not necessary to separate a specific state
        bin_values = limits(i,1):limits(i,2);
    end

    if ~isempty(counts)
        bin_values = find(ismember(counts,bin_values));
    end
    
    % Finally average it 
    average(i) = nanmean(data(bin_values));
end
end

%% Function to sum the counting of something in time bins
function counting = time_bin_counting(data, bin_length, block_length, max_l_expected, n_epochs, state_classification, consider_rest, normalize)

% Check whether is necessary to consider the resting elements (so the last
% bin has a lower number of epochs when compared to the rest)
if isempty(consider_rest)
    consider_rest = false;
end

% Define the beginning and end of each bin
b_length_block = (bin_length*3600)/block_length;
% The last bin considers a lower number of epochs
if consider_rest
    n_bins = ceil(n_epochs/b_length_block);
    beg = ((1:n_bins)-1)*b_length_block + 1;
    en = [(1:n_bins-1)*b_length_block n_epochs];
    limits = [beg' en'];
% The uneven bin is excluded, so all of them have the same length
else
    n_bins = floor(n_epochs/b_length_block);
    beg = ((1:n_bins)-1)*b_length_block + 1;
    en = (1:n_bins)*b_length_block;
    limits = [beg' en'];
end


% If the maximum length has been defined
if ~isempty(max_l_expected)
    if consider_rest
        counting = nan(1,ceil(max_l_expected/b_length_block));
    else
        counting = nan(1,floor(max_l_expected/b_length_block));
    end
else % If there is not an maximum epoch number expected
    counting = nan(1,n_bins);
end


% Define the limits of each bin
for i = 1:n_bins
    % Firstly check whether is necessary to separate a specific state
    if ~isempty(state_classification)
        bin_values = intersect(limits(i,1):limits(i,2),state_classification);
    else % If it is not necessary to separate a specific state
        bin_values = limits(i,1):limits(i,2);
    end
    
    % Finally count it 
    counting(i) = length(find(ismember(data,bin_values)));

    % Normalize by the time in a state
    if normalize && ~isempty(state_classification)
        counting(i) = counting(i)/(length(bin_values) * block_length);
    end
end
end