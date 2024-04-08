% Get coherence for each day

function getspectralcoherence_function(directory,names)
s_p = Signal_Processing;

% Define the band list and band names
band_list = [1 4; 6 9; 4 7; 8 12; 30 80];
band_names = {"Delta",'Theta','Theta2','Theta1','Gamma'};

mkdir(fullfile(directory,'spectral_coherence12'))
mkdir(fullfile(directory,'spectral_coherence13'))
for ii = 1:length(names)
    
    clearvars -except directory names counter ii s_p band_names band_list
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
    
    [Cxy_all12, C_band12, F] = s_p.spectral_coherence(LFP1,LFP2,fs,band_list,band_names);%mudei para LFP3
    sav = fullfile(directory,'spectral_coherence12',char(names{ii}));
    save(sav,"Cxy_all12",'C_band12','F')

    [Cxy_all13, C_band13, F] = s_p.spectral_coherence(LFP1,LFP3,fs,band_list,band_names);%mudei para LFP3
    sav = fullfile(directory,'spectral_coherence13',char(names{ii}));
    save(sav,"Cxy_all13",'C_band13','F')
    

    %[Cxy_all, C_band, F] = s_p.spectral_coherence(LFP1,LFP3,fs,band_list,band_names);%mudei para LFP3
    %sav = fullfile(directory,'spectral_coherence',char(names{ii}));
    %save(sav,"Cxy_all",'C_band','F')
    
    
    disp(names{ii})
    
end



end