%% Complete Routine

% directory = str(base directory name)
directory = 'E:\Barnes Maze - Mestrad\dados matlab\blocked_data';
% names = cell(names of each subdirectory separetelly)
% names = {'B1_D1','B1_D2','B1_D3','B1_D4','B1_D5','B2_D1','B2_D2','B2_D3','B2_D4','B2_D5',...
%     'B3_D1','B3_D2','B3_D3','B3_D4','B3_D5',...
%     'B4_D1','B4_D2','B4_D3','B4_D4','B4_D5','B6_D1','B6_D2','B6_D3','B6_D4','B6_D5'};

names = {'B2_D1','B2_D2','B2_D3','B2_D4','B2_D5',...
    'B3_D1','B3_D2','B3_D3','B3_D4','B3_D5',...
    'B4_D1','B4_D2','B4_D3','B4_D4','B4_D5','B6_D1','B6_D2','B6_D3','B6_D4','B6_D5'};
% names2 = cell(matrix with each subdirectory separatelly
names2 = {'B1_D1','B1_D2','B1_D3','B1_D4','B1_D5';'B2_D1','B2_D2','B2_D3','B2_D4','B2_D5';...
    'B3_D1','B3_D2','B3_D3','B3_D4','B3_D5';...
    'B4_D1','B4_D2','B4_D3','B4_D4','B4_D5';'B6_D1','B6_D2','B6_D3','B6_D4','B6_D5'};

% disp('PSD and Bands')
get_psd_and_band_function(directory,names)
disp('Plot PSD and Bands')
plot_psds_function(directory,names2)
plot_psds_compare_days_function(directory,names2)
plot_bands_throughout_days_function(directory,names,5)
disp('Sleep parameters')
% % FALTANDO
get_sleep_parameters_function(directory,names)
% disp('Spectral coherence')
% getspectralcoherence_function(directory,names)
disp('Phase coherence')
getphasecoherence_function(directory,names)
disp('Plot Spectral and Phase coherence')
% Retirar o plot dos outros 3 dias
% Refazer os plots com base nos diferentes canais de CA1
plot_spectral_phase_coherence_function(directory,names2)
plot_coherence_along_days_function(directory,names,channel_selection)

close all
disp('Ripples')
% Faltando o notchfilter
detect_ripples_block_by_block_function(directory,names)
extract_RDS_info(directory,names)
disp('Delta')
% Faltando detect_Delta_Naty_algorithm_rem
detect_delta_block_by_block_function(directory,names)
detect_delta_block_by_block_CA1_function(directory,names)
disp('Spindle')
% Faltando eegfilt2
detect_spindle_block_by_block_function(directory,names)
disp('Detect RDS coupling')
detect_RDS_block_by_block_function(directory,names,names2)    
disp('Modulation Index')
get_mi_for_each_channel_function(directory, names)
get_mi_for_each_channel_function_n_blocks(directory, names, 3)
get_MI_delta_degree_amplitude(directory,names2,time_limit,10)
get_MI_delta_degree_amplitude_epoch_is_N(directory,names2,time_limit,10)

% disp('Statistics')
channel_selection = [3 3 3 2 3];
time_limit = 180; % Minutes
[results, normalized_results, anova_results, anova_rejected, filtered_ttest, filtered_correlation, results_names] = statistics_function(directory,names2,channel_selection,time_limit,10);
[results, normalized_results, anova_results, anova_rejected, filtered_ttest, filtered_correlation, results_names] = statistics_function_rec_end(directory,names2,channel_selection,time_limit,10);

[results_distribution, anova_results, anova_rejected, results_names, ...
    results_distribution_bin, KW_results, KW_rejected, results_names_bin] = distribution_statistics_function(directory,names2,channel_selection,time_limit,10);

[results_distribution, anova_results, anova_rejected, results_names] = distribution_statistics_function_rec_end(directory,names2,channel_selection,time_limit,10);
[results, normalized_results, anova_results, anova_rejected, correlation_results, filtered_correlation, results_names] = time_bins_statistics_function(directory,names2,channel_selection,[],10,1);
[results, normalized_results, anova_results, anova_rejected, filtered_correlation, results_names] = curve_fit_time_bins_statistics_function(directory,names2,channel_selection,[],10,1);

[curve_fit_results,curve_fit_anova_results,curve_fit_anova_rejected, curve_fit_results_names, ttest_results, filtered_ttest] = run_curve_fit(results, false, 'fit');

% % Create prism files for the new results
dir_to_save = 'E:\Barnes Maze - Mestrad\dados matlab\PrismResults_hour-hour_non_normalized';

create_prism_files_for_results(results_names(136), results, dir_to_save,5)
create_prism_files_for_results(results_names([172;173;176;177]), normalized_results, dir_to_save,4)
create_prism_files_for_results(results_names([174;175;178;179;180]), results, dir_to_save,0)
create_prism_files_for_two_factor_analysis(results_names([134]), results, dir_to_save,5 , 6)
create_prism_files_for_KW_analysis(results_names_bin, results_distribution_bin, dir_to_save,5 , 6)
create_prism_files_for_results_curve_fit(curve_fit_results_names([118 121 122 124 127 134 138 142 119 121 122 123 124 131 146 149 152 155 158]), curve_fit_results, dir_to_save,5)
% First day only
create_prism_files_for_results_first_day(results_names, results, dir_to_save)