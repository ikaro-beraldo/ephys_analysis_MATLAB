% Detect Delta block by block
function detect_delta_block_by_block_CA1_function(directory,names)

mkdir(fullfile(directory,'delta_CA1'))
for ii = 1:length(names)
    disp(names{ii})
    clearvars -except directory names counter ii s_p
    counter = 1;

    % Complete filename
    filename1 = fullfile(directory,names{ii},'blocked_data.mat');
    filename2 = fullfile(directory,names{ii},'GMM_Classification.mat');
    load(filename1,'LFP2','LFP3','fs')
    load(filename2)
    %     LFP1 = gpuArray(LFP1);
    %     LFP2 = gpuArray(LFP2);
    %     LFP3 = gpuArray(LFP3);

    GMM_NREM_All_Sort(:) = false;
    GMM_NREM_All_Sort(GMM.All_Sort == 2) = true;
    GMM_REM_All_Sort(:) = false;
    GMM_REM_All_Sort(GMM.All_Sort == 1) = true;
    GMM_WK_All_Sort(:) = false;
    GMM_WK_All_Sort(GMM.All_Sort == 3) = true;

    % Get only the NREM epochs
    LFP2 = LFP2(GMM.All_Sort == 2,:);
    LFP3 = LFP3(GMM.All_Sort == 2,:);

    % Detect delta
    sd_threshold = 2.8;

    % NATY'S VERSION
    data_windowing.data_w = LFP2;
    data_windowing.time_w = linspace(0,size(LFP2,2)/fs,size(LFP2,2));
    delta_parameters_2 = deltaDetect(data_windowing, fs, sd_threshold, false);
    delta_blocks_2 = [delta_parameters_2.delta_detected(:,1) delta_parameters_2.delta_detected(:,5:7)];
    NREM = find(GMM.All_Sort == 2);
    delta_blocks_2(:,1) = NREM(delta_blocks_2(:,1));

    clear LFP2

    % NATY'S VERSION
    data_windowing.data_w = LFP3;
    data_windowing.time_w = linspace(0,size(LFP3,2)/fs,size(LFP3,2));
    delta_parameters_3 = deltaDetect(data_windowing, fs, sd_threshold, false);
    delta_blocks_3 = [delta_parameters_3.delta_detected(:,1) delta_parameters_3.delta_detected(:,5:7)];
    NREM = find(GMM.All_Sort == 2);
    delta_blocks_3(:,1) = NREM(delta_blocks_3(:,1));

    clear LFP3

    % IKARO'S VERSION
    %     [delta_blocks, parameters] = detect_Delta_Naty_algorithm(LFP1,fs,sd_threshold,false);
    % Fix the delta_blocks epochs according to NREM periods
    %     NREM = find(GMM.All_Sort == 2);
    %     delta_blocks(:,1) = NREM(delta_blocks(:,1));

    display(['total number of candidate events are ' num2str(size(delta_blocks_2,1))])
    display(['total number of candidate events are ' num2str(size(delta_blocks_3,1))])

    sav = fullfile(directory,'delta_CA1',names{ii});
    save(sav,'delta_blocks_2','delta_blocks_3','delta_parameters_2','delta_parameters_3','sd_threshold')

end
end