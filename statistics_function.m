%% Statistics function
function [results, normalized_results, anova_results, anova_rejected, filtered_ttest, filtered_correlation, results_names] = statistics_function(directory,names,channel_selection,time_limit,b_length)

counter_loop = 0;
results = struct;

% Time limit to consider in minutes
time_limit = time_limit * 60 / b_length;

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
        results.bands.LFP1.Delta.NREM(jj,ii) = nanmean(BANDS.Delta.LFP1(NREM));
        results.bands.LFP1.Delta.REM(jj,ii) = nanmean(BANDS.Delta.LFP1(REM));
        results.bands.LFP1.Delta.WK(jj,ii) = nanmean(BANDS.Delta.LFP1(WK));

        % mPFC - Theta
        results.bands.LFP1.Theta.NREM(jj,ii) = nanmean(BANDS.Theta.LFP1(NREM));
        results.bands.LFP1.Theta.REM(jj,ii) = nanmean(BANDS.Theta.LFP1(REM));
        results.bands.LFP1.Theta.WK(jj,ii) = nanmean(BANDS.Theta.LFP1(WK));

        % mPFC - Gamma
        results.bands.LFP1.Gamma.NREM(jj,ii) = nanmean(BANDS.Gamma.LFP1(NREM));
        results.bands.LFP1.Gamma.REM(jj,ii) = nanmean(BANDS.Gamma.LFP1(REM));
        results.bands.LFP1.Gamma.WK(jj,ii) = nanmean(BANDS.Gamma.LFP1(WK));

        % CA11 - Delta
        results.bands.LFP2.Delta.NREM(jj,ii) = nanmean(BANDS.Delta.LFP2(NREM));
        results.bands.LFP2.Delta.REM(jj,ii) = nanmean(BANDS.Delta.LFP2(REM));
        results.bands.LFP2.Delta.WK(jj,ii) = nanmean(BANDS.Delta.LFP2(WK));

        % CA11 - Theta
        results.bands.LFP2.Theta.NREM(jj,ii) = nanmean(BANDS.Theta.LFP2(NREM));
        results.bands.LFP2.Theta.REM(jj,ii) = nanmean(BANDS.Theta.LFP2(REM));
        results.bands.LFP2.Theta.WK(jj,ii) = nanmean(BANDS.Theta.LFP2(WK));

        % CA11 - Gamma
        results.bands.LFP2.Gamma.NREM(jj,ii) = nanmean(BANDS.Gamma.LFP2(NREM));
        results.bands.LFP2.Gamma.REM(jj,ii) = nanmean(BANDS.Gamma.LFP2(REM));
        results.bands.LFP2.Gamma.WK(jj,ii) = nanmean(BANDS.Gamma.LFP2(WK));

        % CA12 - Delta
        results.bands.LFP3.Delta.NREM(jj,ii) = nanmean(BANDS.Delta.LFP3(NREM));
        results.bands.LFP3.Delta.REM(jj,ii) = nanmean(BANDS.Delta.LFP3(REM));
        results.bands.LFP3.Delta.WK(jj,ii) = nanmean(BANDS.Delta.LFP3(WK));

        % CA12 - Theta
        results.bands.LFP3.Theta.NREM(jj,ii) = nanmean(BANDS.Theta.LFP3(NREM));
        results.bands.LFP3.Theta.REM(jj,ii) = nanmean(BANDS.Theta.LFP3(REM));
        results.bands.LFP3.Theta.WK(jj,ii) = nanmean(BANDS.Theta.LFP3(WK));

        % CA12 - Gamma
        results.bands.LFP3.Gamma.NREM(jj,ii) = nanmean(BANDS.Gamma.LFP3(NREM));
        results.bands.LFP3.Gamma.REM(jj,ii) = nanmean(BANDS.Gamma.LFP3(REM));
        results.bands.LFP3.Gamma.WK(jj,ii) = nanmean(BANDS.Gamma.LFP3(WK));

        %% 2 - Script to extract phase coherence throughout the days

        % Load files
        if channel_selection(ii) == 2
            filename = fullfile(directory,'phase_coherence12',char(names{ii,jj}));
        elseif channel_selection(ii) == 3
            filename = fullfile(directory,'phase_coherence13',char(names{ii,jj}));
        end

        load(filename,'Bands')

        % Delta
        results.phase.Delta.NREM(jj,ii) = nanmean(Bands.Delta(NREM));
        results.phase.Delta.REM(jj,ii) = nanmean(Bands.Delta(REM));
        results.phase.Delta.WK(jj,ii) = nanmean(Bands.Delta(WK));

        % Theta
        results.phase.Theta.NREM(jj,ii) = nanmean(Bands.Theta(NREM));
        results.phase.Theta.REM(jj,ii) = nanmean(Bands.Theta(REM));
        results.phase.Theta.WK(jj,ii) = nanmean(Bands.Theta(WK));

        % Gamma
        results.phase.Gamma.NREM(jj,ii) = nanmean(Bands.Gamma(NREM));
        results.phase.Gamma.REM(jj,ii) = nanmean(Bands.Gamma(REM));
        results.phase.Gamma.WK(jj,ii) = nanmean(Bands.Gamma(WK));


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
        results.spectral.Delta.NREM(jj,ii) = nanmean(C_band.Delta(NREM));
        results.spectral.Delta.REM(jj,ii) = nanmean(C_band.Delta(REM));
        results.spectral.Delta.WK(jj,ii) = nanmean(C_band.Delta(WK));

        % Theta
        results.spectral.Theta.NREM(jj,ii) = nanmean(C_band.Theta(NREM));
        results.spectral.Theta.REM(jj,ii) = nanmean(C_band.Theta(REM));
        results.spectral.Theta.WK(jj,ii) = nanmean(C_band.Theta(WK));

        % Gamma
        results.spectral.Gamma.NREM(jj,ii) = nanmean(C_band.Gamma(NREM));
        results.spectral.Gamma.REM(jj,ii) = nanmean(C_band.Gamma(REM));
        results.spectral.Gamma.WK(jj,ii) = nanmean(C_band.Gamma(WK));


        %% 4 - Script to extract MI throughout the days

        % Load files
        filename = fullfile(directory,'MI',char(names{ii,jj}));
        load(filename,'MI_all')

        Pf = [5 12];
        % Amplitude
        Af = [30 55; 65 90; 90 140; 150 170; 180 200; 200 250];

        % LFP2 - 30 - 55
        results.MI.LFP2.f_30_55.NREM(jj,ii) = nanmean(MI_all.LFP2(1,NREM));
        results.MI.LFP2.f_30_55.REM(jj,ii) = nanmean(MI_all.LFP2(1,REM));
        results.MI.LFP2.f_30_55.WK(jj,ii) = nanmean(MI_all.LFP2(1,WK));

        % LFP2 - 65 - 90
        results.MI.LFP2.f_65_90.NREM(jj,ii) = nanmean(MI_all.LFP2(2,NREM));
        results.MI.LFP2.f_65_90.REM(jj,ii) = nanmean(MI_all.LFP2(2,REM));
        results.MI.LFP2.f_65_90.WK(jj,ii) = nanmean(MI_all.LFP2(2,WK));

        % LFP2 - 90 - 140
        results.MI.LFP2.f_90_140.NREM(jj,ii) = nanmean(MI_all.LFP2(3,NREM));
        results.MI.LFP2.f_90_140.REM(jj,ii) = nanmean(MI_all.LFP2(3,REM));
        results.MI.LFP2.f_90_140.WK(jj,ii) = nanmean(MI_all.LFP2(3,WK));

        % LFP2 - 150 - 170
        results.MI.LFP2.f_150_170.NREM(jj,ii) = nanmean(MI_all.LFP2(4,NREM));
        results.MI.LFP2.f_150_170.REM(jj,ii) = nanmean(MI_all.LFP2(4,REM));
        results.MI.LFP2.f_150_170.WK(jj,ii) = nanmean(MI_all.LFP2(4,WK));

        % LFP2 - 180 - 200
        results.MI.LFP2.f_180_200.NREM(jj,ii) = nanmean(MI_all.LFP2(5,NREM));
        results.MI.LFP2.f_180_200.REM(jj,ii) = nanmean(MI_all.LFP2(5,REM));
        results.MI.LFP2.f_180_200.WK(jj,ii) = nanmean(MI_all.LFP2(5,WK));

        % LFP2 - 200 - 250
        results.MI.LFP2.f_200_250.NREM(jj,ii) = nanmean(MI_all.LFP2(6,NREM));
        results.MI.LFP2.f_200_250.REM(jj,ii) = nanmean(MI_all.LFP2(6,REM));
        results.MI.LFP2.f_200_250.WK(jj,ii) = nanmean(MI_all.LFP2(6,WK));


        % LFP3 - 30 - 55
        results.MI.LFP3.f_30_55.NREM(jj,ii) = nanmean(MI_all.LFP3(1,NREM));
        results.MI.LFP3.f_30_55.REM(jj,ii) = nanmean(MI_all.LFP3(1,REM));
        results.MI.LFP3.f_30_55.WK(jj,ii) = nanmean(MI_all.LFP3(1,WK));

        % LFP3 - 65 - 90
        results.MI.LFP3.f_65_90.NREM(jj,ii) = nanmean(MI_all.LFP3(2,NREM));
        results.MI.LFP3.f_65_90.REM(jj,ii) = nanmean(MI_all.LFP3(2,REM));
        results.MI.LFP3.f_65_90.WK(jj,ii) = nanmean(MI_all.LFP3(2,WK));

        % LFP3 - 90 - 140
        results.MI.LFP3.f_90_140.NREM(jj,ii) = nanmean(MI_all.LFP3(3,NREM));
        results.MI.LFP3.f_90_140.REM(jj,ii) = nanmean(MI_all.LFP3(3,REM));
        results.MI.LFP3.f_90_140.WK(jj,ii) = nanmean(MI_all.LFP3(3,WK));

        % LFP3 - 150 - 170
        results.MI.LFP3.f_150_170.NREM(jj,ii) = nanmean(MI_all.LFP3(4,NREM));
        results.MI.LFP3.f_150_170.REM(jj,ii) = nanmean(MI_all.LFP3(4,REM));
        results.MI.LFP3.f_150_170.WK(jj,ii) = nanmean(MI_all.LFP3(4,WK));

        % LFP3 - 180 - 200
        results.MI.LFP3.f_180_200.NREM(jj,ii) = nanmean(MI_all.LFP3(5,NREM));
        results.MI.LFP3.f_180_200.REM(jj,ii) = nanmean(MI_all.LFP3(5,REM));
        results.MI.LFP3.f_180_200.WK(jj,ii) = nanmean(MI_all.LFP3(5,WK));

        % LFP3 - 200 - 250
        results.MI.LFP3.f_200_250.NREM(jj,ii) = nanmean(MI_all.LFP3(6,NREM));
        results.MI.LFP3.f_200_250.REM(jj,ii) = nanmean(MI_all.LFP3(6,REM));
        results.MI.LFP3.f_200_250.WK(jj,ii) = nanmean(MI_all.LFP3(6,WK));

        %% MI considering stable sleep states

        filename = fullfile(directory,'MI_state_block',char(names{ii,jj}));
        load(filename,'MI_all')

        if channel_selection(ii) == 2
            MI_all = MI_all.LFP2;
        else
            MI_all = MI_all.LFP3;
        end

        % LFP3 - 30 - 55
        results.MI_stable.REM.f_30_55(jj,ii) = nanmean(MI_all(1,:));
        %         results.MI_stable.REM.f_30_55.REM(jj,ii) = nanmean(MI_all.LFP3(1,REM));
        %         results.MI_stable.REM.f_30_55.WK(jj,ii) = nanmean(MI_all.LFP3(1,WK));

        % LFP3 - 65 - 90
        results.MI_stable.REM.f_65_90.NREM(jj,ii) = nanmean(MI_all(2,:));
        %         results.MI_stable.REM.f_65_90.REM(jj,ii) = nanmean(MI_all.LFP3(2,REM));
        %         results.MI_stable.REM.f_65_90.WK(jj,ii) = nanmean(MI_all.LFP3(2,WK));

        % LFP3 - 90 - 140
        results.MI_stable.REM.f_90_140.NREM(jj,ii) = nanmean(MI_all(3,:));
        %         results.MI_stable.REM.f_90_140.REM(jj,ii) = nanmean(MI_all.LFP3(3,REM));
        %         results.MI_stable.REM.f_90_140.WK(jj,ii) = nanmean(MI_all.LFP3(3,WK));

        % LFP3 - 150 - 170
        results.MI_stable.REM.f_150_170.NREM(jj,ii) = nanmean(MI_all(4,:));
        %         results.MI_stable.REM.f_150_170.REM(jj,ii) = nanmean(MI_all.LFP3(4,REM));
        %         results.MI_stable.REM.f_150_170.WK(jj,ii) = nanmean(MI_all.LFP3(4,WK));

        % LFP3 - 180 - 200
        results.MI_stable.REM.f_180_200.NREM(jj,ii) = nanmean(MI_all(5,:));
        %         results.MI_stable.REM.f_180_200.REM(jj,ii) = nanmean(MI_all.LFP3(5,REM));
        %         results.MI_stable.REM.f_180_200.WK(jj,ii) = nanmean(MI_all.LFP3(5,WK));

        % LFP3 - 200 - 250
        results.MI_stable.REM.f_200_250.NREM(jj,ii) = nanmean(MI_all(6,:));
        %         results.MI_stable.REM.f_200_250.REM(jj,ii) = nanmean(MI_all.LFP3(6,REM));
        %         results.MI_stable.REM.f_200_250.WK(jj,ii) = nanmean(MI_all.LFP3(6,WK));

        %% Delta MI throught days

        % Load files
        filename = fullfile(directory,'MI_delta',char(names{ii,jj}));
        load(filename,'MI_all')

        Pf = [5 12];
        % Amplitude
        Af = [30 55; 65 90; 90 140; 150 170; 180 200; 200 250];

        % LFP2 - 30 - 55
        results.MI_delta.LFP2.f_30_55.NREM(jj,ii) = nanmean(MI_all.LFP2(1,NREM));
        results.MI_delta.LFP2.f_30_55.REM(jj,ii) = nanmean(MI_all.LFP2(1,REM));
        results.MI_delta.LFP2.f_30_55.WK(jj,ii) = nanmean(MI_all.LFP2(1,WK));

        % LFP2 - 65 - 90
        results.MI_delta.LFP2.f_65_90.NREM(jj,ii) = nanmean(MI_all.LFP2(2,NREM));
        results.MI_delta.LFP2.f_65_90.REM(jj,ii) = nanmean(MI_all.LFP2(2,REM));
        results.MI_delta.LFP2.f_65_90.WK(jj,ii) = nanmean(MI_all.LFP2(2,WK));

        % LFP2 - 90 - 140
        results.MI_delta.LFP2.f_90_140.NREM(jj,ii) = nanmean(MI_all.LFP2(3,NREM));
        results.MI_delta.LFP2.f_90_140.REM(jj,ii) = nanmean(MI_all.LFP2(3,REM));
        results.MI_delta.LFP2.f_90_140.WK(jj,ii) = nanmean(MI_all.LFP2(3,WK));

        % LFP2 - 150 - 170
        results.MI_delta.LFP2.f_150_170.NREM(jj,ii) = nanmean(MI_all.LFP2(4,NREM));
        results.MI_delta.LFP2.f_150_170.REM(jj,ii) = nanmean(MI_all.LFP2(4,REM));
        results.MI_delta.LFP2.f_150_170.WK(jj,ii) = nanmean(MI_all.LFP2(4,WK));

        % LFP2 - 180 - 200
        results.MI_delta.LFP2.f_180_200.NREM(jj,ii) = nanmean(MI_all.LFP2(5,NREM));
        results.MI_delta.LFP2.f_180_200.REM(jj,ii) = nanmean(MI_all.LFP2(5,REM));
        results.MI_delta.LFP2.f_180_200.WK(jj,ii) = nanmean(MI_all.LFP2(5,WK));

        % LFP2 - 200 - 250
        results.MI_delta.LFP2.f_200_250.NREM(jj,ii) = nanmean(MI_all.LFP2(6,NREM));
        results.MI_delta.LFP2.f_200_250.REM(jj,ii) = nanmean(MI_all.LFP2(6,REM));
        results.MI_delta.LFP2.f_200_250.WK(jj,ii) = nanmean(MI_all.LFP2(6,WK));


        % LFP3 - 30 - 55
        results.MI_delta.LFP3.f_30_55.NREM(jj,ii) = nanmean(MI_all.LFP3(1,NREM));
        results.MI_delta.LFP3.f_30_55.REM(jj,ii) = nanmean(MI_all.LFP3(1,REM));
        results.MI_delta.LFP3.f_30_55.WK(jj,ii) = nanmean(MI_all.LFP3(1,WK));

        % LFP3 - 65 - 90
        results.MI_delta.LFP3.f_65_90.NREM(jj,ii) = nanmean(MI_all.LFP3(2,NREM));
        results.MI_delta.LFP3.f_65_90.REM(jj,ii) = nanmean(MI_all.LFP3(2,REM));
        results.MI_delta.LFP3.f_65_90.WK(jj,ii) = nanmean(MI_all.LFP3(2,WK));

        % LFP3 - 90 - 140
        results.MI_delta.LFP3.f_90_140.NREM(jj,ii) = nanmean(MI_all.LFP3(3,NREM));
        results.MI_delta.LFP3.f_90_140.REM(jj,ii) = nanmean(MI_all.LFP3(3,REM));
        results.MI_delta.LFP3.f_90_140.WK(jj,ii) = nanmean(MI_all.LFP3(3,WK));

        % LFP3 - 150 - 170
        results.MI_delta.LFP3.f_150_170.NREM(jj,ii) = nanmean(MI_all.LFP3(4,NREM));
        results.MI_delta.LFP3.f_150_170.REM(jj,ii) = nanmean(MI_all.LFP3(4,REM));
        results.MI_delta.LFP3.f_150_170.WK(jj,ii) = nanmean(MI_all.LFP3(4,WK));

        % LFP3 - 180 - 200
        results.MI_delta.LFP3.f_180_200.NREM(jj,ii) = nanmean(MI_all.LFP3(5,NREM));
        results.MI_delta.LFP3.f_180_200.REM(jj,ii) = nanmean(MI_all.LFP3(5,REM));
        results.MI_delta.LFP3.f_180_200.WK(jj,ii) = nanmean(MI_all.LFP3(5,WK));

        % LFP3 - 200 - 250
        results.MI_delta.LFP3.f_200_250.NREM(jj,ii) = nanmean(MI_all.LFP3(6,NREM));
        results.MI_delta.LFP3.f_200_250.REM(jj,ii) = nanmean(MI_all.LFP3(6,REM));
        results.MI_delta.LFP3.f_200_250.WK(jj,ii) = nanmean(MI_all.LFP3(6,WK));


        %% 5 - Script to extract RDS throughout the days

        % Load files
        filename = fullfile(directory,'ripple',char(names{ii,jj}));
        filename2 = fullfile(directory,'RDS',char(names{ii,jj}));

        load(filename,'ripple_accepted')
        load(filename2,'block')
        ripple_blocks = ripple_accepted(:,1);
        RDS_blocks.delta_ripple = block.delta_ripple_timestamps.delta(:,1);
        RDS_blocks.ripple_delta = block.ripple_delta_timestamps.ripple(:,1);
        RDS_blocks.delta_spindle = block.delta_spindle_timestamps.delta(:,1);
        RDS_blocks.ripple_delta_spindle = block.ripple_delta_spindle_timestamps.ripple(:,1);

        %         if channel_selection(ii) == 2
        %             load(filename,'ripple_result_block_2')
        %             load(filename2,'block_2')
        %             ripple_block = ripple_result_block_2(:,1);
        %             RDS_blocks.delta_ripple = block_2.delta_ripple_timestamps.delta(:,1);
        %             RDS_blocks.ripple_delta = block_2.ripple_delta_timestamps.ripple(:,1);
        %             RDS_blocks.delta_spindle = block_2.delta_spindle_timestamps.delta(:,1);
        %             RDS_blocks.ripple_delta_spindle = block_2.ripple_delta_spindle_timestamps.ripple(:,1);
        %         elseif channel_selection(ii) == 3
        %             load(filename,'ripple_result_block_3')
        %             load(filename2,'block_3')
        %             ripple_blocks = ripple_result_block_3(:,1);
        %             RDS_blocks.delta_ripple = block_3.delta_ripple_timestamps.delta(:,1);
        %             RDS_blocks.ripple_delta = block_3.ripple_delta_timestamps.ripple(:,1);
        %             RDS_blocks.delta_spindle = block_3.delta_spindle_timestamps.delta(:,1);
        %             RDS_blocks.ripple_delta_spindle = block_3.ripple_delta_spindle_timestamps.ripple(:,1);
        %         end

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

        results.RDS.ripple(jj,ii) = 60*length(find(ripple_blocks <= time_limit)) / (length(NREM)*b_length);
        results.RDS.delta(jj,ii) = 60*length(find(delta_blocks <= time_limit)) / (length(NREM)*b_length);
        results.RDS.spindle(jj,ii) = 60*length(find(spindle_blocks <= time_limit)) / (length(NREM)*b_length);
        results.RDS.ripple_delta(jj,ii) = 60*length(find(RDS_blocks.ripple_delta <= time_limit)) / (length(NREM)*b_length);
        results.RDS.delta_ripple(jj,ii) = 60*length(find(RDS_blocks.delta_ripple <= time_limit)) / (length(NREM)*b_length);
        results.RDS.delta_spindle(jj,ii) = 60*length(find(RDS_blocks.delta_spindle <= time_limit)) / (length(NREM)*b_length);
        results.RDS.ripple_delta_spindle(jj,ii) = 60*length(find(RDS_blocks.ripple_delta_spindle <= time_limit)) / (length(NREM)*b_length);
        results.RDS.delta_CA11(jj,ii) = 60*length(find(delta_blocks_2 <= time_limit)) / (length(NREM)*b_length);
        results.RDS.delta_CA12(jj,ii) = 60*length(find(delta_blocks_3 <= time_limit)) / (length(NREM)*b_length);

        %% Extract RDS parameters

        % RIPPLE
        filename = fullfile(directory,'ripple',char(names{ii,jj}));
        load(filename,'ripple_result_block_2','ripple_result_block_3','duration','RMS','PEAK','FREQ')
        time_select_2 = find(ripple_result_block_2(:,1) <= time_limit);
        time_select_3 = find(ripple_result_block_3(:,1) <= time_limit);
        time_select = find(ripple_accepted(:,1) <= time_limit);

        results.RDS.duration.ripple_2(jj,ii) = mean(duration.ripple_2(time_select_2));
        results.RDS.duration.ripple_3(jj,ii) = mean(duration.ripple_3(time_select_3));
        results.RDS.RMS.ripple_2(jj,ii) = mean(RMS.ripple_2(time_select_2));
        results.RDS.RMS.ripple_3(jj,ii) = mean(RMS.ripple_3(time_select_3));
        results.RDS.PEAK.ripple_2(jj,ii) = mean(PEAK.ripple_2(time_select_2));
        results.RDS.PEAK.ripple_3(jj,ii) = mean(PEAK.ripple_3(time_select_3));
        results.RDS.FREQ.ripple_2(jj,ii) = mean(FREQ.ripple_2(time_select_2));
        results.RDS.FREQ.ripple_3(jj,ii) = mean(FREQ.ripple_3(time_select_3));

        % Get a unique ripple (based on selected channels)
        results.RDS.duration.ripple(jj,ii) = mean(duration.ripple(time_select));
        dur_ripple = duration.ripple(time_select);
        results.RDS.duration.ripple_counts_higher_than_100ms(jj,ii) = length(find(dur_ripple>=0.1))/length(dur_ripple);
        results.RDS.RMS.ripple(jj,ii) = mean(RMS.ripple(time_select));
        results.RDS.PEAK.ripple(jj,ii) = mean(PEAK.ripple(time_select));
        results.RDS.FREQ.ripple(jj,ii) = mean(FREQ.ripple(time_select));

%         % Get a unique ripple (based on selected channels)
%         if channel_selection(ii) == 2
%             results.RDS.duration.ripple(jj,ii) = mean(duration.ripple_2(time_select_2));
%             results.RDS.RMS.ripple(jj,ii) = mean(RMS.ripple_2(time_select_2));
%             results.RDS.PEAK.ripple(jj,ii) = mean(PEAK.ripple_2(time_select_2));
%             results.RDS.FREQ.ripple(jj,ii) = mean(FREQ.ripple_2(time_select_2));
%         elseif channel_selection(ii) == 3
%             results.RDS.duration.ripple(jj,ii) = mean(duration.ripple_3(time_select_3));
%             results.RDS.RMS.ripple(jj,ii) = mean(RMS.ripple_3(time_select_3));
%             results.RDS.PEAK.ripple(jj,ii) = mean(PEAK.ripple_3(time_select_3));
%             results.RDS.FREQ.ripple(jj,ii) = mean(FREQ.ripple_3(time_select_3));
%         end


        % SPINDLE
        filename = fullfile(directory,'spindle',char(names{ii,jj}));
        load(filename,'spindle_blocks','duration','RMS','PEAK','FREQ')
        time_select_2 = find(spindle_blocks(:,1) <= time_limit);

        results.RDS.duration.spindle(jj,ii) = mean(duration.spindle(time_select_2));
        results.RDS.RMS.spindle(jj,ii) = mean(RMS.spindle(time_select_2));
        results.RDS.PEAK.spindle(jj,ii) = mean(PEAK.spindle(time_select_2));
        results.RDS.FREQ.spindle(jj,ii) = mean(FREQ.spindle(time_select_2));

        % DELTA
        filename = fullfile(directory,'delta',char(names{ii,jj}));
        load(filename,'delta_blocks','delta_parameters')
        time_select_2 = find(delta_blocks(:,1) <= time_limit);

        results.RDS.Amplitude_onset.delta(jj,ii) = mean(delta_parameters.parameters(time_select_2,2));
        results.RDS.duration.delta(jj,ii) = mean(delta_parameters.parameters(time_select_2,3));
        results.RDS.Slope_onset.delta(jj,ii) = mean(delta_parameters.parameters(time_select_2,4));
        results.RDS.Slope_offset.delta(jj,ii) = mean(delta_parameters.parameters(time_select_2,5));
        results.RDS.Area.delta(jj,ii) = mean(delta_parameters.parameters(time_select_2,6));
        results.RDS.Energy.delta(jj,ii) = mean(delta_parameters.parameters(time_select_2,7));

        % DELTA CA1
        filename = fullfile(directory,'delta_CA1',char(names{ii,jj}));
        load(filename,'delta_blocks_2','delta_parameters_2','delta_blocks_3','delta_parameters_3')
        time_select_2 = find(delta_blocks_2(:,1) <= time_limit);
        time_select_3 = find(delta_blocks_3(:,1) <= time_limit);

        % Delta CA11
        results.RDS.Amplitude_onset.delta_CA11(jj,ii) = mean(delta_parameters_2.parameters(time_select_2,2));
        results.RDS.duration.delta_CA11(jj,ii) = mean(delta_parameters_2.parameters(time_select_2,3));
        results.RDS.Slope_onset.delta_CA11(jj,ii) = mean(delta_parameters_2.parameters(time_select_2,4));
        results.RDS.Slope_offset.delta_CA11(jj,ii) = mean(delta_parameters_2.parameters(time_select_2,5));
        results.RDS.Area.delta_CA11(jj,ii) = mean(delta_parameters_2.parameters(time_select_2,6));
        results.RDS.Energy.delta_CA11(jj,ii) = mean(delta_parameters_2.parameters(time_select_2,7));

        % Delta CA12
        results.RDS.Amplitude_onset.delta_CA12(jj,ii) = mean(delta_parameters_3.parameters(time_select_3,2));
        results.RDS.duration.delta_CA12(jj,ii) = mean(delta_parameters_3.parameters(time_select_3,3));
        results.RDS.Slope_onset.delta_CA12(jj,ii) = mean(delta_parameters_3.parameters(time_select_3,4));
        results.RDS.Slope_offset.delta_CA12(jj,ii) = mean(delta_parameters_3.parameters(time_select_3,5));
        results.RDS.Area.delta_CA12(jj,ii) = mean(delta_parameters_3.parameters(time_select_3,6));
        results.RDS.Energy.delta_CA12(jj,ii) = mean(delta_parameters_3.parameters(time_select_3,7));

        

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


%% Normalize data by the first day

% Get every single results field
results_names = fieldnamesr(results,'prefix');

for ii = 1:length(results_names)
    cm = sprintf('matrix = %s;',results_names{ii});
    % execute command
    eval(cm)

    % Normalize for the first day
    matrix = matrix./matrix(1,:);
    cm = sprintf('normalized_%s = matrix;', results_names{ii});
    eval(cm)

    % Normalize zscore
    %     matrix = zscore(matrix,[],1);
    %     cm = sprintf('normalized_%s = matrix;', results_names{ii});
    %     eval(cm)
end

%% Calculate anova

% Get every single results field
results_names = fieldnamesr(results,'prefix');
% Create a table for anova results
anova_results = cell(length(results_names),1);

for ii = 1:length(results_names)

    % execute command
    eval(sprintf("[anova_results{ii,1},anova_results{ii,2},anova_results{ii,3}] = anova1(%s',[],'off');",results_names{ii}))

end

% Get the results names
anova_results = [anova_results results_names];

% Get the comparisons which wielded a p-value < 0.05
rejected = find([anova_results{:,1}] < 0.05);
anova_rejected = {results_names{rejected}; anova_results{rejected,1}; anova_results{rejected,2}; anova_results{rejected,3}}';

%% Calculate ttest2

% Get every single results field
results_names = fieldnamesr(normalized_results,'prefix');
% Create a table for anova results
ttest2_results = cell(length(results_names),1);

for ii = 1:length(results_names)
    eval(sprintf("matrix = %s';",results_names{ii}));
    for jj = 1:size(matrix,2)
        for kk = 1:size(matrix,2)
            eval(sprintf("a = %s';",results_names{ii}))
            a = a(:,jj);
            eval(sprintf("b = %s';",results_names{ii}))
            b = b(:,kk);
            % execute command
            [results_tt_h(jj,kk), results_tt_p(jj,kk)] = ttest(a,b);
            %eval(sprintf("[ttest2_results{ii,1},ttest2_results{ii,2}] = ttest2(%s',[],'off');",results_names{ii}))
        end
    end
    ttest_results{ii,1} = results_tt_h;
    ttest_results{ii,2} = results_tt_p;
    ttest_results{ii,3} = results_names{ii};
    clear results_tt_h results_tt_p
end

% Filter results
filtered_ttest = {};
count = 1;
for ii = 1:size(ttest_results,1)
    if nansum(find(ttest_results{ii,1} == 1)) > 1
        filtered_ttest{count,1} = ttest_results{ii,1};
        filtered_ttest{count,2} = ttest_results{ii,2};
        filtered_ttest{count,3} = ttest_results{ii,3};
        count = count + 1;
    end
end

% Get the comparisons which wielded a p-value < 0.05
% rejected = find([anova_results{:,1}] < 0.05);
% anova_rejected = {results_names{rejected} anova_results(rejected,1) anova_results(rejected,2) anova_results(rejected,3)};


%% CORRELATION PROB TRIAL STATISTICS
% Load the Barnes Maze Prob Trial statistics
load("E:\Barnes Maze - Mestrad\Resultados DLC - Barnes\ProbTrial\Final_results.mat")
prob_time_on_each_quadrant = task_parameters.("Time_each_quadrant");

% Get every single results field
results_names = fieldnamesr(normalized_results,'prefix');
correlation_results_prob = cell(length(results_names),4);
counter = 1;
for jj = 1:length(results_names)
    % Get the parameters
    a = prob_time_on_each_quadrant(:,1)./120;
    eval(sprintf('b = %s;',results_names{jj}));

    % Take the test parameter and average its days, except the baseline
    b = nanmean(b(2:end,:));

    [correlation_results_prob{counter,1},correlation_results_prob{counter,2}] = corrcoef(a,b);
    % Fix
    correlation_results_prob{counter,1} = correlation_results_prob{counter,1}(2);
    correlation_results_prob{counter,2} = correlation_results_prob{counter,2}(2);
    correlation_results_prob{counter,3} = 'prob_time_each_quadrant';
    correlation_results_prob{counter,4} = results_names{jj};

    counter = counter + 1;
end

%% Calculate Pearson Correlation for each combination

% Get every single results field
results_names = fieldnamesr(results,'prefix');
correlation_results = cell(length(results_names)*length(results_names),4);
counter = 1;
for ii = 1:length(results_names)
    for jj = 1:length(results_names)
        eval(sprintf('a = %s;',results_names{ii}));
        eval(sprintf('b = %s;',results_names{jj}));

        % Check if a and b have the same length
        if sum(size(a) == size(b)) < 2
            if size(a,1) > size(b,1)
                a = a(2:5,:);
            else
                b = b(2:5,:);
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