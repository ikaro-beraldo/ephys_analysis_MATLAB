%%time_resolved_phase_amplitude_coupling
clearvars -except id
clc

name = 'B1_D1';

%Select the files saved from Parkinson Sleep script
% disp('Select the saved files from ParkinsonSleep script')
% text = sprintf('Select the saved files from ParkinsonSleep script');
% [baseName, folder] = uigetfile('',text);
% file = fullfile(folder, baseName);
% %Load selected file
% load(file,'LFP2_block','n_All_Sort','FS','B')

load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data',name,'blocked_data.mat'),'LFP3','fs')
load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data',name,'GMM_Classification.mat'),'GMM')

% if isempty(id)
%     id = input('Insert the animal id (ex: C3A3_29_11_2018): ','s');
%     id = string(id);
% end

lfp2 = LFP3;   %HIPPO
n_All_Sort = GMM.All_Sort;
B = 10;
block_length = 10;
FS = fs;

clear LFP3

%Transpose data to avoid errors
% if length(lfp2) == FS*B
%     lfp2 = transpose(lfp2);
% end

data_length = length(lfp2(1,:));
data_heigth_total = size(lfp2,1);   %######## Defines the number of periods being analyzed,
%also defining the time spent
srate = FS;

%%

clearvars -except n_All_Sort PhaseFreqVector AmpFreqVector ComodulogramMeanHH ...
    ComodulogramMeanHP ComodulogramHH ComodulogramHP ...
    data_length lfp1 lfp2 srate id data_heigth_total folder block_length

n_All_Sort(data_heigth_total+1:end) = 0;

REM = find(n_All_Sort == 1);

%STATES = vector with the period numbers for each states
STATES = REM;
data_heigth = length(STATES);

AmpFreqVector=20:10:250;

PhaseFreq_BandWidth = 2;
AmpFreq_BandWidth= 20;
% AmpFreqVector=20:10:310;

%% For comodulation calculation (only has to be calculated once)
nbin = 18;
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for c=1:nbin
    position(c) = -pi+(c-1)*winsize;
end
%% Do filtering and Hilbert transform on CPU

% Start and end periods
start_period = 1;
end_period =  10; %Standard --> length(STATES)

ComodulogramHH=single(zeros(length(AmpFreqVector),block_length-3,length(start_period:end_period)));
selected_freq = zeros(length(AmpFreqVector),block_length-3,length(start_period:end_period));
theta_power = zeros(29,block_length-3,length(start_period:end_period));
% ComodulogramHP=single(zeros(length(1),length(AmpFreqVector),data_heigth));
% MeanAmpTotal1 = zeros(1,nbin,data_heigth);
% MeanAmpTotal2 = zeros(1,nbin,data_heigth);

%disp('CPU filtering')
fprintf('Periods: 1 to %d: ',data_heigth_total);
for blck = start_period:end_period%length(STATES)
    
    fprintf('%d . ', blck);
    
    %% Define the Amplitude- and Phase- Frequencies
    
    %The phase frequency is defined by the theta frequency which peaks
    %(bandwidth = 2hz)
    
    %     %Filtered signal in gamma components
    %     gamma_filtered = eegfilt2(lfp2(STATES(blck),:),srate,20,310);
    %     %Hilbert transform
    %     gamma_filtered = abs(hilbert(gamma_filtered));
    %     %Find the power peak frequency of the envelope of the signal
    %     L = length(gamma_filtered);
    %     Y = abs(fft(gamma_filtered)/L);
    %     Y = Y(1:L/2+1);
    %     Y(2:end-1) = 2*Y(2:end-1);
    %     freq = srate*(0:(L/2))/L;
    %
    %     %Ignore frequency components under 0.5hz
    %     accept = find(freq > 0.5);
    %
    %
    %     [~,freq_peak] = max(Y(accept));
    %
    %     PhaseFreqVector = freq(accept(freq_peak));
    
    %%
    
    AmpFreqTransformedHH = zeros(length(AmpFreqVector), data_length);
    
    for ii=1:length(AmpFreqVector)
        Af1 = AmpFreqVector(ii);
        Af2=Af1+AmpFreq_BandWidth;
        AmpFreqHH = eegfilt2(lfp2(STATES(blck),:),srate,Af1,Af2); % just filtering
        AmpFreqTransformedHH = abs(hilbert(AmpFreqHH)); % getting the amplitude envelope
        
        %Segmentation of data (4 sec with 75% overlap)
        sample_idx = 1:block_length*srate;
        segmented_idx = buffer(sample_idx,4*srate,4*srate*0.75);
        segmented_idx(:,find(segmented_idx(1,:) == 0)) = [];
        
        for seg = 1:size(segmented_idx,2)
            %             a = (seg-1)*srate+1;
            %             b = (seg + 3)*srate;
            signal = AmpFreqTransformedHH(segmented_idx(:,seg));
            signal = signal - mean(signal);
                
           
            %Find the power peak frequency of the envelope of the signal
            L = length(signal);
            Y = abs(fft(signal)/L);
            Y = Y(1:L/2+1);
            Y(2:end-1) = 2*Y(2:end-1);
            freq = srate*(0:(L/2))/L;
            
            %Get theta power in each 4 sec block
            if ii == 1
                [pxx,f] = pwelch(lfp2(STATES(blck),segmented_idx(:,seg)),srate,0.5,srate*4,srate);
                theta = find(f >=5 & f <=12);
                theta_power(:,seg,blck) = pxx(theta);
            end
            %             Y(find(freq < 0.5)) = 0;
            
            %Ignore frequency components under 0.5hz
            %     accept = find(freq > 0.5);
            %     [~,freq_peak] = max(Y(accept));
            
            limits = find(freq >= 5 & freq <= 9);
            [~,freq_peak] = max(Y(limits));
            freq_peak = freq(limits(freq_peak));
            
            %Selected phase frequency
            selected_freq(ii,seg,blck) = freq_peak;
            
            PhaseFreqVector = freq_peak;
            Pf1 = PhaseFreqVector;
            Pf2 = PhaseFreqVector + PhaseFreq_BandWidth;
%             Pf1 = PhaseFreqVector - PhaseFreq_BandWidth/2;
%             Pf2 = PhaseFreqVector + PhaseFreq_BandWidth/2;
%             PhaseFreq = eegfilt2(lfp2(STATES(blck),:),srate,Pf1,Pf2);
%             PhaseFreq = PhaseFreq(segmented_idx(:,seg));
            PhaseFreq=eegfilt2(lfp2(STATES(blck),segmented_idx(:,seg)),srate,Pf1,Pf2); % this is just filtering
            PhaseFreqTransformed = angle(hilbert(PhaseFreq));
            
            [MI1,~]=ModIndex_v2(PhaseFreqTransformed, signal, position);
            ComodulogramHH(ii,seg,blck)=MI1;
            
        end
    end
    
    %     for jj=1:length(PhaseFreqVector)
    %         Pf1 = PhaseFreqVector(jj);
    %         Pf2 = Pf1 + PhaseFreq_BandWidth;
    %         PhaseFreq=eegfilt2(lfp2(STATES(blck),:),srate,Pf1,Pf2); % this is just filtering
    %         PhaseFreqTransformed(jj,:) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
    %     end
end

%% Plotting

%Getting the maximum values of each period
for period = 1:blck
    max_periods(period) = abs(max(ComodulogramHH(:,:,period),[],'all'));
    min_periods(period) = abs(min(ComodulogramHH(:,:,period),[],'all'));
end

%Plot each period
for period = 1:blck
    subplot(3,1,1)
    pcolor(1:block_length-3,AmpFreqVector,abs(ComodulogramHH(:,:,period)))
     %surf(1:block_length,AmpFreqVector,abs(ComodulogramHH(:,:,period)),'EdgeColor','none')
    % axis xy; axis tight; colormap(jet); view(0,90);
     %contourf(1:block_length-3,AmpFreqVector,abs(ComodulogramHH(:,:,period)),100,'lines','none');
    % s0.FaceColor = 'interp';
%     xlabel('Time (s)')
    ylabel('Amplitude Frequency (Hz)')
    h = colorbar;
    ylabel(h, 'MI')
    caxis([mean(min_periods) min(max_periods)])
    colormap(jet)
    title('Phase-Amplitude Coupling')
    
    subplot(3,1,2)
    contourf(1:block_length-3,AmpFreqVector,abs(selected_freq(:,:,period)),100,'lines','none');
    % s1.FaceColor = 'interp';
%     xlabel('Time (s)')
    ylabel('Amplitude Frequency (Hz)')
    h = colorbar;
    ylabel(h, 'Phase Frequency (Hz)')
    title('Selected Phase Frequency')
    colormap(jet)
    
    subplot(3,1,3)
%     contourf(1:block_length-3,f(theta),theta_power(:,:,period),100,'lines','none');
    pcolor(1:block_length-3,f(theta),theta_power(:,:,period));
shading interp
    % s0.FaceColor = 'interp';
    xlabel('Time (s)')
    ylabel('Theta Frequency (Hz)')
    h = colorbar;
    ylabel(h, 'Power (mV²)') 
    title('Theta Power')
    % caxis([mean(min_periods) mean(max_periods)])
    colormap(jet)
    pause

end

%Plot all the periods concatenated
figure
comodulogram_concatenated = reshape(ComodulogramHH,size(ComodulogramHH,1),size(ComodulogramHH,2)*size(ComodulogramHH,3),1);
pcolor(1:size(comodulogram_concatenated,2),AmpFreqVector,abs(comodulogram_concatenated));
shading flat 
% s.FaceColor = 'interp';
% s.LineWidth = 0.1;
colormap(jet)
colorbar
caxis([0 min(max_periods)])
ylabel('Amplitude Frequency (Hz)')
xlabel('Time (s)')
h = colorbar;
ylabel(h, 'MI')
title('REM periods')

