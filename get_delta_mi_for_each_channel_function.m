%% Modulation Index

% Get coherence for each day
function get_delta_mi_for_each_channel_function(directory, names)

s_p = Signal_Processing;

mkdir(fullfile(directory,'MI_delta'))
for ii = 1:length(names)
    fprintf('%s\n.',names{ii})
    clearvars -except directory names counter ii s_p
    counter = 1;
    
    % Complete filename
    filename1 = fullfile(directory,names{ii},'blocked_data.mat');
    filename2 = fullfile(directory,names{ii},'GMM_Classification.mat');
    load(filename1)
    load(filename2)
    % Fix All_Sort
    GMM_NREM_All_Sort(:) = false;
    GMM_NREM_All_Sort(GMM.All_Sort == 2) = true;
    GMM_REM_All_Sort(:) = false;
    GMM_REM_All_Sort(GMM.All_Sort == 1) = true;
    GMM_WK_All_Sort(:) = false;
    GMM_WK_All_Sort(GMM.All_Sort == 3) = true;
    
    % Theta
    Pf = [0.9 4];
    % Amplitude
    Af = [30 55; 65 90; 90 140; 150 170; 180 200; 200 250];
    
    % Pre-allocate variables
    MI_all.LFP2 = nan(size(Af,1),size(LFP2,1));
    MI_all.LFP3 = nan(size(Af,1),size(LFP3,1));
    
    MeanAmp_all.LFP2 = nan(size(Af,1),size(LFP2,1),18);
    MeanAmp_all.LFP3 = nan(size(Af,1),size(LFP3,1),18);
    
    % Modulation Index calculation
    for ep = 1:size(LFP3,1)
        fprintf('%d.',ep)
        for bb = 1:size(Af,1)
            % Calculate MI for each amplitude band
            % LFP3
            [MI,MeanAmp] = ModIndex_v1(LFP2(ep,:),fs,Pf(1),Pf(2),Af(bb,1),Af(bb,2));
            MI_all.LFP2(bb,ep) = MI;
            MeanAmp_all.LFP2(bb,ep,:) = MeanAmp;
            
            % LFP2
            [MI,MeanAmp] = ModIndex_v1(LFP3(ep,:),fs,Pf(1),Pf(2),Af(bb,1),Af(bb,2));
            MI_all.LFP3(bb,ep) = MI;
            MeanAmp_all.LFP3(bb,ep,:) = MeanAmp;
        end
    end
    
    sav = fullfile(directory,'MI_delta',char(names{ii}));
    save(sav,'MeanAmp_all','MI_all','Pf','Af')

    fprintf('\n')
end
end

