% Get coherence for each day

function getphasecoherence_function(directory,names)
s_p = Signal_Processing;

mkdir(fullfile(directory,'phase_coherence12'))
mkdir(fullfile(directory,'phase_coherence13'))
for ii = 1:length(names)
    
    mkdir(fullfile(directory,'PSD'))
    clearvars -except directory names counter ii s_p
    counter = 1;
    
    % Complete filename
    filename1 = fullfile(directory,names{ii},'blocked_data.mat');
    filename2 = fullfile(directory,names{ii},'GMM_Classification.mat');
    load(filename1)
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
    
    [Call12, phiall12, confCall12,f] = s_p.phase_coherence(LFP1,LFP2,fs);%%mudei para LFP3
    sav = fullfile(directory,'phase_coherence12',char(names{ii}));
    save(sav,"Call12",'phiall12','confCall12','f')

    [Call13, phiall13, confCall13,f] = s_p.phase_coherence(LFP1,LFP3,fs);%%mudei para LFP3
    sav = fullfile(directory,'phase_coherence13',char(names{ii}));
    save(sav,"Call13",'phiall13','confCall13','f')

%     [Call, phiall, confCall,f] = s_p.phase_coherence(LFP1,LFP3,fs);%%mudei para LFP3
%     sav = fullfile(directory,'phase_coherence',char(names{ii}));
%     save(sav,"Call",'phiall','confCall','f')
    
    
    
    disp(names{ii})
    
end

%% Get the phase coherence for each defined band (LFP1 LFP2)

s_p = Signal_Processing;

for ii = 1:length(names)
    
    clearvars -except directory names counter ii s_p
    counter = 1;
    
    % Complete filename
    filename1 = fullfile(directory,'phase_coherence12',names{ii});
    
    load(filename1,'Call12','f');
    
    % Separate bands
    delta = find(f >=0.9 & f<=4);
    theta = find(f >= 6 & f<=9);
    gamma = find(f > 30);
    
    % Get the average;
    Bands.Delta = nanmean(Call12(delta,:));
    Bands.Theta = nanmean(Call12(theta,:));
    Bands.Gamma = nanmean(Call12(gamma,:));
    
    sav = fullfile(directory,'phase_coherence12',char(names{ii}));
    save(sav,'Bands','-append')
    
    
end

%% %% Get the phase coherence for each defined band (LFP1 LFP3)

s_p = Signal_Processing;

for ii = 1:length(names)
    
    clearvars -except directory names counter ii s_p
    counter = 1;
    
    % Complete filename
    filename1 = fullfile(directory,'phase_coherence13',names{ii});
    
    load(filename1,'Call13','f');
    
    % Separate bands
    delta = find(f >=0.9 & f<=4);
    theta = find(f >= 6 & f<=9);
    gamma = find(f > 30);
    
    % Get the average;
    Bands.Delta = nanmean(Call13(delta,:));
    Bands.Theta = nanmean(Call13(theta,:));
    Bands.Gamma = nanmean(Call13(gamma,:));
    
    sav = fullfile(directory,'phase_coherence13',char(names{ii}));
    save(sav,'Bands','-append')
    
    
end
end
