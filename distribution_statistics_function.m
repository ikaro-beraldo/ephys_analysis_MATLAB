function [results_distribution, anova_results, anova_rejected, results_names,...
    results_distribution_bin, KW_results, KW_rejected, results_names_bin]...
    = distribution_statistics_function(directory,names,channel_selection,time_limit,b_length)

counter_loop = 0;
results_distribution = struct;
%% Result distribution
% Pre-allocate
results_distribution.RDS.duration.ripple_2 = cell(size(names,2),1);
results_distribution.RDS.duration.ripple_3 = cell(size(names,2),1);
results_distribution.RDS.RMS.ripple_2 = cell(size(names,2),1);
results_distribution.RDS.RMS.ripple_3 = cell(size(names,2),1);
results_distribution.RDS.PEAK.ripple_2 = cell(size(names,2),1);
results_distribution.RDS.PEAK.ripple_3 = cell(size(names,2),1);
results_distribution.RDS.FREQ.ripple_2 = cell(size(names,2),1);
results_distribution.RDS.FREQ.ripple_3 = cell(size(names,2),1);

results_distribution.RDS.duration.ripple = cell(size(names,2),1);
results_distribution.RDS.RMS.ripple = cell(size(names,2),1);
results_distribution.RDS.PEAK.ripple = cell(size(names,2),1);
results_distribution.RDS.FREQ.ripple = cell(size(names,2),1);

% Spectral and Phase coherecence
pot = {'Delta','Theta','Gamma'};
state = {'NREM','REM','WK'};
analysis = {'phase','spectral'};
for x = 1:length(analysis)
    for y = 1:length(pot)
        for z = 1:length(state)
            eval(sprintf('results_distribution.%s.%s.%s = cell(size(names,2),1);',analysis{x},pot{y},state{z}))
        end
    end
end

%% Result distribution BIN
n_bins = 7;
hour_limit = 6;

% Pre-allocate
results_distribution_bin.RDS.duration.ripple_2 = cell(size(names,2),n_bins);
results_distribution_bin.RDS.duration.ripple_3 = cell(size(names,2),n_bins);
results_distribution_bin.RDS.RMS.ripple_2 = cell(size(names,2),n_bins);
results_distribution_bin.RDS.RMS.ripple_3 = cell(size(names,2),n_bins);
results_distribution_bin.RDS.PEAK.ripple_2 = cell(size(names,2),n_bins);
results_distribution_bin.RDS.PEAK.ripple_3 = cell(size(names,2),n_bins);
results_distribution_bin.RDS.FREQ.ripple_2 = cell(size(names,2),n_bins);
results_distribution_bin.RDS.FREQ.ripple_3 = cell(size(names,2),n_bins);

results_distribution_bin.RDS.duration.ripple = cell(size(names,2),n_bins);
results_distribution_bin.RDS.RMS.ripple = cell(size(names,2),n_bins);
results_distribution_bin.RDS.PEAK.ripple = cell(size(names,2),n_bins);
results_distribution_bin.RDS.FREQ.ripple = cell(size(names,2),n_bins);

% Spectral and Phase coherecence
pot = {'Delta','Theta','Gamma'};
state = {'NREM','REM','WK'};
analysis = {'phase','spectral'};
for x = 1:length(analysis)
    for y = 1:length(pot)
        for z = 1:length(state)
            eval(sprintf('results_distribution_bin.%s.%s.%s = cell(size(names,2),n_bins);',analysis{x},pot{y},state{z}))
        end
    end
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


%%
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

        % For Bin analysis
        NREM_bin = NREM;
        REM_bin = REM;
        WK_bin = WK;

        % Filter by time
        NREM = NREM(NREM <= time_limit);
        REM = REM(REM <= time_limit);
        WK = WK(WK <= time_limit);


        %% Ripple features by day - ignoring animal identities

        % RIPPLE
        filename = fullfile(directory,'ripple',char(names{ii,jj}));
        load(filename,'ripple_result_block_2','ripple_result_block_3','ripple_accepted','duration','RMS','PEAK','FREQ')
        time_select_2 = find(ripple_result_block_2(:,1) <= time_limit);
        time_select_3 = find(ripple_result_block_3(:,1) <= time_limit);
        time_select = find(ripple_accepted(:,1) <= time_limit);

        % Get the results
        results_distribution.RDS.duration.ripple_2{jj,1} = [results_distribution.RDS.duration.ripple_2{jj,1}; duration.ripple_2(time_select_2)];
        results_distribution.RDS.duration.ripple_3{jj,1} = [results_distribution.RDS.duration.ripple_3{jj,1}; duration.ripple_3(time_select_3)];
        results_distribution.RDS.RMS.ripple_2{jj,1} = [results_distribution.RDS.RMS.ripple_2{jj,1}; RMS.ripple_2(time_select_2)'];
        results_distribution.RDS.RMS.ripple_3{jj,1} = [results_distribution.RDS.RMS.ripple_3{jj,1}; RMS.ripple_3(time_select_3)'];
        results_distribution.RDS.PEAK.ripple_2{jj,1} = [results_distribution.RDS.PEAK.ripple_2{jj,1}; PEAK.ripple_2(time_select_2)'];
        results_distribution.RDS.PEAK.ripple_3{jj,1} = [results_distribution.RDS.PEAK.ripple_3{jj,1}; PEAK.ripple_3(time_select_3)'];
        results_distribution.RDS.FREQ.ripple_2{jj,1} = [results_distribution.RDS.FREQ.ripple_2{jj,1}; FREQ.ripple_2(time_select_2)'];
        results_distribution.RDS.FREQ.ripple_3{jj,1} = [results_distribution.RDS.FREQ.ripple_3{jj,1}; FREQ.ripple_3(time_select_3)'];

        results_distribution.RDS.duration.ripple{jj,1} = [results_distribution.RDS.duration.ripple{jj,1}; duration.ripple(time_select)];
        results_distribution.RDS.RMS.ripple{jj,1} = [results_distribution.RDS.RMS.ripple{jj,1}; RMS.ripple(time_select)'];
        results_distribution.RDS.PEAK.ripple{jj,1} = [results_distribution.RDS.PEAK.ripple{jj,1}; PEAK.ripple(time_select)'];
        results_distribution.RDS.FREQ.ripple{jj,1} = [results_distribution.RDS.FREQ.ripple{jj,1}; FREQ.ripple(time_select)'];

        % Get the data distribution and concatenate cells element-wise
        get_data = time_bin_get(duration.ripple_2, 1, 10, max_rec_length, NREM_bin, true, ripple_result_block_2(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.duration.ripple_2(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.duration.ripple_2(jj,:), get_data, 'UniformOutput', false);
        get_data = time_bin_get(duration.ripple_3, 1, 10, max_rec_length, NREM_bin, true, ripple_result_block_3(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.duration.ripple_3(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.duration.ripple_3(jj,:), get_data, 'UniformOutput', false);

        get_data = time_bin_get(RMS.ripple_2, 1, 10, max_rec_length, NREM_bin, true, ripple_result_block_2(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.RMS.ripple_2(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.RMS.ripple_2(jj,:), get_data, 'UniformOutput', false);
        get_data = time_bin_get(RMS.ripple_3, 1, 10, max_rec_length, NREM_bin, true, ripple_result_block_3(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.RMS.ripple_3(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.RMS.ripple_3(jj,:), get_data, 'UniformOutput', false);

        get_data = time_bin_get(PEAK.ripple_2, 1, 10, max_rec_length, NREM_bin, true, ripple_result_block_2(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.PEAK.ripple_2(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.PEAK.ripple_2(jj,:), get_data, 'UniformOutput', false);
        get_data = time_bin_get(PEAK.ripple_3, 1, 10, max_rec_length, NREM_bin, true, ripple_result_block_3(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.PEAK.ripple_3(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.PEAK.ripple_3(jj,:), get_data, 'UniformOutput', false);

        get_data = time_bin_get(FREQ.ripple_2, 1, 10, max_rec_length, NREM_bin, true, ripple_result_block_2(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.FREQ.ripple_2(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.FREQ.ripple_2(jj,:), get_data, 'UniformOutput', false);
        get_data = time_bin_get(FREQ.ripple_3, 1, 10, max_rec_length, NREM_bin, true, ripple_result_block_3(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.FREQ.ripple_3(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.FREQ.ripple_3(jj,:), get_data, 'UniformOutput', false);

        % Get the data distribution and concatenate cells element-wise
        get_data = time_bin_get(duration.ripple, 1, 10, max_rec_length, NREM_bin, true, ripple_accepted(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.duration.ripple(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.duration.ripple(jj,:), get_data, 'UniformOutput', false);
        
        get_data = time_bin_get(RMS.ripple, 1, 10, max_rec_length, NREM_bin, true, ripple_accepted(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.RMS.ripple(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.RMS.ripple(jj,:), get_data, 'UniformOutput', false);
        
        get_data = time_bin_get(PEAK.ripple, 1, 10, max_rec_length, NREM_bin, true, ripple_accepted(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.PEAK.ripple(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.PEAK.ripple(jj,:), get_data, 'UniformOutput', false);
       
        get_data = time_bin_get(FREQ.ripple, 1, 10, max_rec_length, NREM_bin, true, ripple_accepted(:,1), length(GMM.All_Sort));
        results_distribution_bin.RDS.FREQ.ripple(jj,:) = cellfun(@vertcat, results_distribution_bin.RDS.FREQ.ripple(jj,:), get_data, 'UniformOutput', false);
      

        %% 2 - Script to extract phase coherence throughout the days

        % Load files
        if channel_selection(ii) == 2
            filename = fullfile(directory,'phase_coherence12',char(names{ii,jj}));
        elseif channel_selection(ii) == 3
            filename = fullfile(directory,'phase_coherence13',char(names{ii,jj}));
        end

        load(filename,'Bands')

        % Delta
        results_distribution.phase.Delta.NREM{jj,1} = [results_distribution.phase.Delta.NREM{jj,1}; Bands.Delta(NREM)'];
        results_distribution.phase.Delta.REM{jj,1} = [results_distribution.phase.Delta.REM{jj,1}; Bands.Delta(REM)'];
        results_distribution.phase.Delta.WK{jj,1} = [results_distribution.phase.Delta.WK{jj,1}; Bands.Delta(WK)'];

        % Theta
        results_distribution.phase.Theta.NREM{jj,1} = [results_distribution.phase.Theta.NREM{jj,1}; Bands.Theta(NREM)'];
        results_distribution.phase.Theta.REM{jj,1} = [results_distribution.phase.Theta.REM{jj,1}; Bands.Theta(REM)'];
        results_distribution.phase.Theta.WK{jj,1} = [results_distribution.phase.Theta.WK{jj,1}; Bands.Theta(WK)'];

        % Gamma
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

        % Delta
        results_distribution.spectral.Delta.NREM{jj,1} = [results_distribution.spectral.Delta.NREM{jj,1}; C_band.Delta(NREM)'];
        results_distribution.spectral.Delta.REM{jj,1} = [results_distribution.spectral.Delta.REM{jj,1}; C_band.Delta(REM)'];
        results_distribution.spectral.Delta.WK{jj,1} = [results_distribution.spectral.Delta.WK{jj,1}; C_band.Delta(WK)'];

        % Theta
        results_distribution.spectral.Theta.NREM{jj,1} = [results_distribution.spectral.Theta.NREM{jj,1}; C_band.Theta(NREM)'];
        results_distribution.spectral.Theta.REM{jj,1} = [results_distribution.spectral.Theta.REM{jj,1}; C_band.Theta(REM)'];
        results_distribution.spectral.Theta.WK{jj,1} = [results_distribution.spectral.Theta.WK{jj,1}; C_band.Theta(WK)'];

        % Gamma
        results_distribution.spectral.Gamma.NREM{jj,1} = [results_distribution.spectral.Gamma.NREM{jj,1}; C_band.Gamma(NREM)'];
        results_distribution.spectral.Gamma.REM{jj,1} = [results_distribution.spectral.Gamma.REM{jj,1}; C_band.Gamma(REM)'];
        results_distribution.spectral.Gamma.WK{jj,1} = [results_distribution.spectral.Gamma.WK{jj,1}; C_band.Gamma(WK)'];


       %% 3.1 - Spectral and Phase coherence - Hour/hour

       % Spectral and Phase coherecence
       pot = {'Delta','Theta','Gamma'};
       state = {'NREM','REM','WK'};
       for y = 1:length(pot)
           for z = 1:length(state)
               % Phase
               eval(sprintf('get_data = time_bin_get(Bands.%s, 1, 10, max_rec_length, %s_bin, true, [], length(GMM.All_Sort));',pot{y},state{z}))
               eval(sprintf('results_distribution_bin.phase.%s.%s(jj,:) = cellfun(@vertcat, results_distribution_bin.phase.%s.%s(jj,:), get_data, "UniformOutput", false);',pot{y},state{z},pot{y},state{z}))
               % Spectral
               eval(sprintf('get_data = time_bin_get(C_band.%s, 1, 10, max_rec_length, %s_bin, true, [], length(GMM.All_Sort));',pot{y},state{z}))
               eval(sprintf('results_distribution_bin.spectral.%s.%s(jj,:) = cellfun(@vertcat, results_distribution_bin.spectral.%s.%s(jj,:), get_data, "UniformOutput", false);',pot{y},state{z},pot{y},state{z}))
           end
       end

      

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


%% Organize results to Kruskall-Wallis

% Get every single results field
results_names_bin = fieldnamesr(results_distribution_bin,'prefix');

% Loop for each result
for ii = 1:length(results_names_bin)
    % Create a cell with group indexing
    groups = cell(size(results_distribution_bin.RDS.duration.ripple_2,1),hour_limit);
    % Pre-allocate the organized_bin_parameter struct
    eval(sprintf('organized_%s = [];',results_names_bin{ii}))
    eval(sprintf('group_%s = [];',results_names_bin{ii}))
    % Get a vector with group indexing for kruskall-wallis
    for jj = 1:numel(groups)
        command1 = sprintf('groups{jj} = jj*ones(length(%s{jj}),1);',results_names_bin{ii},results_names_bin{ii});
        eval(command1)
        % Concatenate all the elements to a single column vector
        command2 = sprintf('organized_%s = [organized_%s; %s{jj}];',results_names_bin{ii},results_names_bin{ii},results_names_bin{ii});
        eval(command2)
        % Concate the group indexing as well
        eval(sprintf('group_%s = [group_%s; groups{jj}];',results_names_bin{ii},results_names_bin{ii}))
    end
end

%% Run Kruskall-Wallis
% Get every single results field
results_names = fieldnamesr(organized_results_distribution_bin,'prefix');
group_values = fieldnamesr(group_results_distribution_bin,'prefix');
% Create a table for anova results
KW_results = cell(length(results_names),1);

for ii = 1:length(results_names)
    % execute command
    eval(sprintf("[KW_results{ii,1},KW_results{ii,2},KW_results{ii,3}] = kruskalwallis(%s,%s,'off');",results_names{ii},group_values{ii}))
end

% Get the results names
KW_results = [KW_results results_names];

% Get the comparisons which wielded a p-value < 0.05
rejected = find([KW_results{:,1}] < 0.05);
KW_rejected = {results_names{rejected}; KW_results{rejected,1}; KW_results{rejected,2}; KW_results{rejected,3}}';


end

%% Function to get the temporal series in time bins
function get_data = time_bin_get(data, bin_length, block_length, max_l_expected, state_classification, consider_rest, counts, n_epochs)

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
        get_data = cell(1,ceil(max_l_expected/b_length_block));
    else
        get_data = cell(1,floor(max_l_expected/b_length_block));
    end
else % If there is not an maximum epoch number expected
    get_data = cell(1,n_bins);
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
    
    % Finally get the data regarding the temporal series time bin 
    % Check whether it is a column vector (IMPORTANT)
    if iscolumn(data(bin_values))
        get_data{i} = data(bin_values);
    else
        get_data{i} = data(bin_values)';
    end
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