%% Modulation Index

% Get coherence for each day
function get_mi_for_each_channel_function_n_blocks(directory, names, n_blocks, only_rem)
% 
% % Run only for REM periods
% if isempty(only_rem)
%     only_rem = false;
% end

s_p = Signal_Processing;

mkdir(fullfile(directory,'MI'))
for ii = 1:length(names)
    fprintf('%s\n',names{ii})
    clearvars -except directory names counter ii s_p n_blocks
    counter = 1;
    
    % Complete filename
    filename1 = fullfile(directory,names{ii},'blocked_data.mat');
    filename2 = fullfile(directory,names{ii},'GMM_Classification.mat');
    load(filename1,'LFP2','LFP3','fs')
    load(filename2)
    % Fix All_Sort
    GMM_NREM_All_Sort(:) = false;
    GMM_NREM_All_Sort(GMM.All_Sort == 2) = true;
    GMM_REM_All_Sort(:) = false;
    GMM_REM_All_Sort(GMM.All_Sort == 1) = true;
    GMM_WK_All_Sort(:) = false;
    GMM_WK_All_Sort(GMM.All_Sort == 3) = true;

    % Define the idx of n_blocks (default = 3)
    [REM_idx, NREM_idx, WK_idx] = define_beg_end_for_n_blocks(GMM.All_Sort,n_blocks);
    
    % Theta
    Pf = [5 9];
    % Amplitude
    Af = [30 55; 65 90; 90 140; 150 170; 180 200; 200 250];
    
    % Pre-allocate variables
    MI_all.LFP2 = nan(size(Af,1),size(REM_idx,1));
    MI_all.LFP3 = nan(size(Af,1),size(REM_idx,1));
    
    MeanAmp_all.LFP2 = nan(size(Af,1),size(REM_idx,1),18);
    MeanAmp_all.LFP3 = nan(size(Af,1),size(REM_idx,1),18);
    
    fprintf('Total: %d\n',size(REM_idx,1))
    % Modulation Index calculation
    for ep = 1:size(REM_idx,1)
        fprintf('%d.',ep)
        for bb = 1:size(Af,1)
            % Define the data block
            a = REM_idx(ep,1);
            b = REM_idx(ep,2);
            data_LFP2 = reshape(LFP2(a:b,:)',[],1)';
            data_LFP3 = reshape(LFP3(a:b,:)',[],1)';
            % Calculate MI for each amplitude band
            % LFP2
            [MI,MeanAmp] = ModIndex_v1(data_LFP2,fs,Pf(1),Pf(2),Af(bb,1),Af(bb,2));
            MI_all.LFP2(bb,ep) = MI;
            MeanAmp_all.LFP2(bb,ep,:) = MeanAmp;
            
            % LFP3
            [MI,MeanAmp] = ModIndex_v1(data_LFP3,fs,Pf(1),Pf(2),Af(bb,1),Af(bb,2));
            MI_all.LFP3(bb,ep) = MI;
            MeanAmp_all.LFP3(bb,ep,:) = MeanAmp;
        end
    end
    
    sav = fullfile(directory,'MI_state_block',char(names{ii}));
    save(sav,'MeanAmp_all','MI_all','Pf','Af','REM_idx')

    fprintf('\n')
end
end

function [REM_idx, NREM_idx, WK_idx] = define_beg_end_for_n_blocks(All_Sort,n_blocks)

    % Get REM blocks
    REM = zeros(1,length(All_Sort));
    REM(All_Sort == 1) = 1;
    REM_idx = get_b_and_g(REM, n_blocks);

    % Get NREM blocks
    NREM = zeros(1,length(All_Sort));
    NREM(All_Sort == 2) = 1;
    NREM_idx = get_b_and_g(NREM, n_blocks);

    % Get WK blocks
    WK = zeros(1,length(All_Sort));
    WK(All_Sort == 1) = 1;
    WK_idx = get_b_and_g(WK, n_blocks);
    
end

% Define the beginning and end of each block
function state_idx = get_b_and_g(state, n_blocks)
primary_detection_b = find(diff(state)>0)+1;
primary_detection_e = find(diff(state)<0);

if primary_detection_b(1)>primary_detection_e(1)
    primary_detection_e = primary_detection_e(2:end);
end

if primary_detection_b(end)>primary_detection_e(end)
    primary_detection_b = primary_detection_b(1:end-1);
end

selected_idx = find(((primary_detection_e - primary_detection_b)+1) >= n_blocks);
state_idx = [primary_detection_b(selected_idx)' primary_detection_e(selected_idx)'];
end

