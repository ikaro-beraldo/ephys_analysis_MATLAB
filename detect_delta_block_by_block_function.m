% Detect Delta block by block
function detect_delta_block_by_block_function(directory,names)

mkdir(fullfile(directory,'delta'))
for ii = 1:length(names)
    disp(names{ii})
    clearvars -except directory names counter ii s_p
    counter = 1;

    % Complete filename
    filename1 = fullfile(directory,names{ii},'blocked_data.mat');
    filename2 = fullfile(directory,names{ii},'GMM_Classification.mat');
    load(filename1,'LFP1','fs')
    load(filename2,'GMM')

    % Filter data
    filtered_LFP1 = eegfilt2(LFP1,fs,1,4);

    % Find the artifact epochs
    [row, ~] = find(filtered_LFP1 >= 1500);
    row = unique(row);

    NREM = find(GMM.All_Sort == 2);
    % Extract the artifact epochs
    [~,artifact] = intersect(NREM,row);
    NREM(artifact) = [];

    filtered_LFP1_aux = filtered_LFP1(NREM,:);

    LFP1_aux = LFP1(NREM,:);

    % Detect delta
    sd_threshold = 2.7;

    % NATY'S VERSION
    data_windowing.data_w = LFP1_aux;
    data_windowing.time_w = linspace(0,size(LFP1_aux,2)/fs,size(LFP1_aux,2));
    [delta_parameters, ~] = deltaDetect_Ikaro(data_windowing, filtered_LFP1_aux, fs, sd_threshold, false);
    delta_blocks = [delta_parameters.delta_detected(:,1) delta_parameters.delta_detected(:,5:7)];
    delta_blocks(:,1) = NREM(delta_blocks(:,1));


    display(['total number of candidate events are ' num2str(size(delta_blocks,1))])

    sav = fullfile(directory,'delta',names{ii});
    save(sav,'delta_blocks','delta_parameters','sd_threshold')

end
end