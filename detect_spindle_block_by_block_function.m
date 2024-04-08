% Detect spindles block by block
function detect_spindle_block_by_block_function(directory,names)

mkdir(fullfile(directory,'spindle'))
for ii = 1:length(names)
    disp(names{ii})
    clearvars -except directory names counter ii s_p
    counter = 1;
    
    % Complete filename
    filename1 = fullfile(directory,names{ii},'blocked_data.mat');
    filename2 = fullfile(directory,names{ii},'GMM_Classification.mat');
    load(filename1,'LFP1','fs')
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
    
    
    % Detect spindle
    spindle_blocks = detectSpindle(LFP1,[],[],GMM.All_Sort,size(LFP1,2)/fs,fs);
    
    %display(['total number of candidate events are ' num2str(size(spindles_blocks,1))])
    
    sav = fullfile(directory,'spindle',names{ii});
    save(sav,'spindle_blocks')
    
    
end
end
