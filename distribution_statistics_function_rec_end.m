function [results_distribution, anova_results, anova_rejected, results_names] = distribution_statistics_function_rec_end(directory,names,channel_selection,time_limit,b_length)

counter_loop = 0;
results_distribution = struct;

% Pre-allocate
results_distribution.RDS.duration.ripple_2 = cell(size(names,2),1);
results_distribution.RDS.duration.ripple_3 = cell(size(names,2),1);
results_distribution.RDS.RMS.ripple_2 = cell(size(names,2),1);
results_distribution.RDS.RMS.ripple_3 = cell(size(names,2),1);
results_distribution.RDS.PEAK.ripple_2 = cell(size(names,2),1);
results_distribution.RDS.PEAK.ripple_3 = cell(size(names,2),1);
results_distribution.RDS.FREQ.ripple_2 = cell(size(names,2),1);
results_distribution.RDS.FREQ.ripple_3 = cell(size(names,2),1);

% Coherence PHASE
results_distribution.phase.Delta.NREM = cell(size(names,2),1);
results_distribution.phase.Delta.REM = cell(size(names,2),1);
results_distribution.phase.Delta.WK = cell(size(names,2),1);
results_distribution.phase.Theta.NREM = cell(size(names,2),1);
results_distribution.phase.Theta.REM = cell(size(names,2),1);
results_distribution.phase.Theta.WK = cell(size(names,2),1);
results_distribution.phase.Gamma.NREM = cell(size(names,2),1);
results_distribution.phase.Gamma.REM = cell(size(names,2),1);
results_distribution.phase.Gamma.WK = cell(size(names,2),1);

% Coherence SPECTRAL
results_distribution.spectral.Delta.NREM = cell(size(names,2),1);
results_distribution.spectral.Delta.REM = cell(size(names,2),1);
results_distribution.spectral.Delta.WK = cell(size(names,2),1);
results_distribution.spectral.Theta.NREM = cell(size(names,2),1);
results_distribution.spectral.Theta.REM = cell(size(names,2),1);
results_distribution.spectral.Theta.WK = cell(size(names,2),1);
results_distribution.spectral.Gamma.NREM = cell(size(names,2),1);
results_distribution.spectral.Gamma.REM = cell(size(names,2),1);
results_distribution.spectral.Gamma.WK = cell(size(names,2),1);

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
        NREM = NREM(NREM >= time_limit);
        REM = REM(REM >= time_limit);
        WK = WK(WK >= time_limit);

        %% Ripple features by day - ignoring animal identities

        % RIPPLE
        filename = fullfile(directory,'ripple',char(names{ii,jj}));
        load(filename,'ripple_result_block_2','ripple_result_block_3','duration','RMS','PEAK','FREQ')
        time_select_2 = find(ripple_result_block_2(:,1) >= time_limit);
        time_select_3 = find(ripple_result_block_3(:,1) >= time_limit);

        % Get the results
        results_distribution.RDS.duration.ripple_2{jj,1} = [results_distribution.RDS.duration.ripple_2{jj,1}; duration.ripple_2(time_select_2)];
        results_distribution.RDS.duration.ripple_3{jj,1} = [results_distribution.RDS.duration.ripple_3{jj,1}; duration.ripple_3(time_select_3)];
        results_distribution.RDS.RMS.ripple_2{jj,1} = [results_distribution.RDS.RMS.ripple_2{jj,1}; RMS.ripple_2(time_select_2)'];
        results_distribution.RDS.RMS.ripple_3{jj,1} = [results_distribution.RDS.RMS.ripple_3{jj,1}; RMS.ripple_3(time_select_3)'];
        results_distribution.RDS.PEAK.ripple_2{jj,1} = [results_distribution.RDS.PEAK.ripple_2{jj,1}; PEAK.ripple_2(time_select_2)];
        results_distribution.RDS.PEAK.ripple_3{jj,1} = [results_distribution.RDS.PEAK.ripple_3{jj,1}; PEAK.ripple_3(time_select_3)];
        results_distribution.RDS.FREQ.ripple_2{jj,1} = [results_distribution.RDS.FREQ.ripple_2{jj,1}; FREQ.ripple_2(time_select_2)'];
        results_distribution.RDS.FREQ.ripple_3{jj,1} = [results_distribution.RDS.FREQ.ripple_3{jj,1}; FREQ.ripple_3(time_select_3)'];

        %% 2 - Script to extract phase coherence throughout the days

        % Load files
        if channel_selection(ii) == 2
            filename = fullfile(directory,'phase_coherence12',char(names{ii,jj}));
        elseif channel_selection(ii) == 3
            filename = fullfile(directory,'phase_coherence13',char(names{ii,jj}));
        end

        load(filename,'Bands')

        % Coherence PHASE
        results_distribution.phase.Delta.NREM{jj,1} = [results_distribution.phase.Delta.NREM{jj,1}; Bands.Delta(NREM)'];
        results_distribution.phase.Delta.REM{jj,1} = [results_distribution.phase.Delta.REM{jj,1}; Bands.Delta(REM)'];
        results_distribution.phase.Delta.WK{jj,1} = [results_distribution.phase.Delta.WK{jj,1}; Bands.Delta(WK)'];
        results_distribution.phase.Theta.NREM{jj,1} = [results_distribution.phase.Theta.NREM{jj,1}; Bands.Theta(NREM)'];
        results_distribution.phase.Theta.REM{jj,1} = [results_distribution.phase.Theta.REM{jj,1}; Bands.Theta(REM)'];
        results_distribution.phase.Theta.WK{jj,1} = [results_distribution.phase.Theta.WK{jj,1}; Bands.Theta(WK)'];
        results_distribution.phase.Gamma.NREM{jj,1} = [results_distribution.phase.Gamma.NREM{jj,1}; Bands.Gamma(NREM)'];
        results_distribution.phase.Gamma.REM{jj,1} = [results_distribution.phase.Gamma.REM{jj,1}; Bands.Gamma(REM)'];
        results_distribution.phase.Gamma.WK{jj,1} = [results_distribution.phase.Gamma.WK{jj,1}; Bands.Gamma(WK)'];


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
        
        % Coherence SPECTRAL
        results_distribution.spectral.Delta.NREM{jj,1} = [results_distribution.spectral.Delta.NREM{jj,1}; C_band.Delta(NREM)'];
        results_distribution.spectral.Delta.REM{jj,1} = [results_distribution.spectral.Delta.REM{jj,1}; C_band.Delta(REM)'];
        results_distribution.spectral.Delta.WK{jj,1} = [results_distribution.spectral.Delta.WK{jj,1}; C_band.Delta(WK)'];
        results_distribution.spectral.Theta.NREM{jj,1} = [results_distribution.spectral.Theta.NREM{jj,1}; C_band.Theta(NREM)'];
        results_distribution.spectral.Theta.REM{jj,1} = [results_distribution.spectral.Theta.REM{jj,1}; C_band.Theta(REM)'];
        results_distribution.spectral.Theta.WK{jj,1} = [results_distribution.spectral.Theta.WK{jj,1}; C_band.Theta(WK)'];
        results_distribution.spectral.Gamma.NREM{jj,1} = [results_distribution.spectral.Gamma.NREM{jj,1}; C_band.Gamma(NREM)'];
        results_distribution.spectral.Gamma.REM{jj,1} = [results_distribution.spectral.Gamma.REM{jj,1}; C_band.Gamma(REM)'];
        results_distribution.spectral.Gamma.WK{jj,1} = [results_distribution.spectral.Gamma.WK{jj,1}; C_band.Gamma(WK)'];


    end
end
%% Organize results to anova

% Get every single results field
results_names = fieldnamesr(results_distribution,'prefix');
% Loop for each result
for ii = 1:length(results_names)
    for jj = 1:size(names,2)
        % Get a vector with group indexing for anova
        command1 = sprintf('%s{jj,2} = jj*ones(length(%s{jj,1}),1);',results_names{ii},results_names{ii});
        eval(command1)
        command2 = sprintf('organized_%s = vertcat(%s{:,1});',results_names{ii},results_names{ii});
        eval(command2)
        command3 = sprintf('group_%s = vertcat(%s{:,2});',results_names{ii},results_names{ii});
        eval(command3)
    end
end

%% Run anova
% Get every single results field
results_names = fieldnamesr(organized_results_distribution,'prefix');
group_values = fieldnamesr(group_results_distribution,'prefix');
% Create a table for anova results
anova_results = cell(length(results_names),1);

for ii = 1:length(results_names)

    % execute command
    eval(sprintf("[anova_results{ii,1},anova_results{ii,2},anova_results{ii,3}] = anova1(%s,%s,'off');",results_names{ii},group_values{ii}))

end

% Get the results names
anova_results = [anova_results results_names];

% Get the comparisons which wielded a p-value < 0.05
rejected = find([anova_results{:,1}] < 0.05);
anova_rejected = {results_names{rejected}; anova_results{rejected,1}; anova_results{rejected,2}; anova_results{rejected,3}}';

end