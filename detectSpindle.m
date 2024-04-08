%Spindle Detector
% MODIFICAÇÕES 24/03 --> smooth de com uma janela de 100 ms e o sinal
% filtrado é elevado ao quadrado antes  

function timestamps = detectSpindle(pfc_data,filtered_data_10_16,filtered_data_20_30,estados,block_length,fs)


%Noise threshold
thresh2 = 0.3;
thresh1 = 2;
%Time threshold [min max] in seconds
time_threshold = [0.4 3];

%Filt the PFC data, case it is not already filtered
if isempty(filtered_data_10_16) == 1
    %If the data is a column vector, transposes it
    if size(pfc_data,1) > size(pfc_data,2)
        pfc_data = pfc_data';
    end
    
    %Bandpass filt the pfc data (10-16hz spindle band)
    filtered_data_10_16 = eegfilt2(pfc_data,fs,10,[]);
    filtered_data_10_16 = eegfilt2(filtered_data_10_16,fs,[],16);

end

%Get the instataneous amplitude by Hilbert transform
%Filt
% fl1 = 300;
% % [hil_filt_10_16,~] =
%
% [hil_filt_10_16,aux] = envelope(filtered_data_10_16,fl1,'analytic');
% clear aux

for v = 1:size(filtered_data_10_16,1)
    hil_filt_10_16(v,:) = abs(smooth(hilbert(filtered_data_10_16(v,:).^2),fs*100));
end

SWS = find(estados == 2);
%Get the SWS samples
SWS_samples = [];

for i = 1:length(SWS)
    a = (SWS(i)-1) * block_length * fs + 1;
    b = SWS(i) * block_length * fs;
    SWS_samples = [SWS_samples a:b];
end

SWS(SWS > size(pfc_data,1)) = [];

%Get the threshold of spindles based only on SWS periods
filt_mean = mean(nanmean(hil_filt_10_16(SWS,:)));
filt_std = std2(hil_filt_10_16(~isnan(hil_filt_10_16)));

%Peak value to be considered an event (mean + 1 * std)
filt_threshold = filt_mean + thresh1 * filt_std;
%Defines the beginning and end of the event (mean + 0.2 * std)
filt_threshold_2 = filt_mean + thresh2 * filt_std;

%Detection loop
trigger = 0;
event = [];
v = 1;
%Block loop
for e = 1:size(hil_filt_10_16,1)
    %Sample loop
    i = 1;
    fprintf('.%d',e);
    
    if e == size(hil_filt_10_16,1)
        fprintf('\n');
    end
    
    while i <= size(hil_filt_10_16,2)
        %Find putative spindles if the envelope is higher than 'peak
        %threshold'
        if hil_filt_10_16(e,i) >= filt_threshold
            %Event matrix [block beginning peak end]. Beginning and end
            %will be added later
            event(v,:) = [e 0 i 0];
            
            z = i-1;
            %Loop to find the beginning
            while z > 0
                %Find the first sample before the peak to be equal to mean
                % + 0.2 Standard deviation
                if hil_filt_10_16(e,z) <= filt_threshold_2
                    event(v,2) = z;
                    z = 0;
                else
                    z = z - 1;
                end
            end
            
            %Indicates that exist an end to the event
            f_end = false;
            
            z = i + 1;
            %Loop to find the end
            while z <= size(hil_filt_10_16,2)
                %Find the first sample before the peak to be equal to mean
                % + 0.2 Standard deviation
                if hil_filt_10_16(e,z) <= filt_threshold_2
                    event(v,4) = z;
                    z = size(hil_filt_10_16,2)+1;
                    jump_event= event(v,4);
                    f_end = true;
                else
                    z = z + 1;
                end
            end
            
            %Increase 1 row to the event matrix
            v = v + 1;
            %Jumps to sample right after the end of an event
            if f_end == true
                i = jump_event;
            end
            
        end
        i = i + 1;
    end
end

%Excluding by time
%Analyse each putative event and find their duration
%Must be between 0.5 and 2.5 seconds
keep_event = [];
v = 1;
for i = 1:size(event,1)
    if event(i,4) - event(i,2) >= time_threshold(1) * fs && event(i,4) - event(i,2) <= time_threshold(2) * fs
        keep_event(v,:) = event(i,:);
        v = v + 1;
    end
end

%Select only SWS periods
keep_sws = [];
for v = 1:length(SWS)
    keep_sws = [keep_sws; find(keep_event(:,1) == SWS(v))];
end

keep_event = keep_event(keep_sws,:);

%Exclude periods with no beginning or end
exc = [find(keep_event(:,2) == 0 ); find(keep_event(:,4) == 0)];
exc = unique(exc);
keep_event(exc,:) = [];

%% Exlude events based on high instantaneous amplitude on 20-30 hz

%Bandpass filt the pfc data(20-30hz subthreshold)
    filtered_data_20_30 = eegfilt2(pfc_data,fs,20,[]);
    filtered_data_20_30 = eegfilt2(filtered_data_20_30,fs,[],30);

%Mean amplitude value of 20-30hz filtered data
per_mean = mean(nanmean(filtered_data_20_30));
% Amplitude value standard deviation of 20-30hz filtered data
per_std = std2(filtered_data_20_30(~isnan(filtered_data_20_30)));

%Threshold (mean + 4.5 * std)
check_thresh = per_mean + 4.5 * per_std;
%Check loop
for v = 1:length(keep_event)
    
    h_value = find(filtered_data_20_30(keep_event(v,1),keep_event(v,2):keep_event(v,4)) > check_thresh);
    
    if isempty(h_value) == 0
        keep_event(v,:) = 0;
    end
    clear h_value
end
%Clear the events in which the instantaneou amplitude values are higher
%than threshold
keep_event(find(keep_event(:,1) == 0),:) = [];

%Total events
fprintf('\nNumber of events: %d\n',length(keep_event));

clear filtered_data_20_30

%% Exlude events based on high instantaneous amplitude on 4-9 hz

 %Bandpass filt the pfc data(4-9hz subthreshold)
    filtered_data_4_9 = eegfilt2(pfc_data,fs,4,[]);
    filtered_data_4_9 = eegfilt2(filtered_data_4_9,fs,[],9);

%Mean amplitude value of 20-30hz filtered data
per_mean = mean(nanmean(filtered_data_4_9));
% Amplitude value standard deviation of 20-30hz filtered data
per_std = std2(filtered_data_4_9(~isnan(filtered_data_4_9)));

%Threshold (mean + 4.5 * std)
check_thresh = per_mean + 4.5 * per_std;
%Check loop
for v = 1:size(keep_event,1)
    
    h_value = find(filtered_data_4_9(keep_event(v,1),keep_event(v,2):keep_event(v,4)) > check_thresh);
    
    if isempty(h_value) == 0
        keep_event(v,:) = 0;
    end
    clear h_value
end
%Clear the events in which the instantaneou amplitude values are higher
%than threshold
keep_event(find(keep_event(:,1) == 0),:) = [];

%Total events
fprintf('\nNumber of events: %d\n',length(keep_event));

clear filtered_data_4_9

%% Exlude events based on high instantaneous amplitude on 17-32 hz

%Bandpass filt the pfc data(17-32 hz)
    filtered_data_17_32 = eegfilt2(pfc_data,fs,17,[]);
    filtered_data_17_32 = eegfilt2(filtered_data_17_32,fs,[],32);

%Mean amplitude value of 20-30hz filtered data
per_mean = mean(nanmean(filtered_data_17_32));
% Amplitude value standard deviation of 20-30hz filtered data
per_std = std2(filtered_data_17_32(~isnan(filtered_data_17_32)));

%Threshold (mean + 4.5 * std)
check_thresh = per_mean + 4.5 * per_std;
%Check loop
for v = 1:size(keep_event,1)
    
    h_value = find(filtered_data_17_32(keep_event(v,1),keep_event(v,2):keep_event(v,4)) > check_thresh);
    
    if isempty(h_value) == 0
        keep_event(v,:) = 0;
    end
    clear h_value
end
%Clear the events in which the instantaneou amplitude values are higher
%than threshold
keep_event(find(keep_event(:,1) == 0),:) = [];

%Total events
fprintf('\nNumber of events: %d\n',length(keep_event));

clear filtered_data_17_32

%% Exlude events based on high instantaneous amplitude on 35_65 hz

 %Bandpass filt the pfc data(35-65 hz)
    filtered_data_35_65 = eegfilt2(pfc_data,fs,35,[]);
    filtered_data_35_65 = eegfilt2(filtered_data_35_65,fs,[],65);

%Mean amplitude value of 20-30hz filtered data
per_mean = mean(nanmean(filtered_data_35_65));
% Amplitude value standard deviation of 20-30hz filtered data
per_std = std2(filtered_data_35_65(~isnan(filtered_data_35_65)));

%Threshold (mean + 4.5 * std)
check_thresh = per_mean + 4.5 * per_std;
%Check loop
for v = 1:size(keep_event,1)
    
    h_value = find(filtered_data_35_65(keep_event(v,1),keep_event(v,2):keep_event(v,4)) > check_thresh);
    
    if isempty(h_value) == 0
        keep_event(v,:) = 0;
    end
    clear h_value
end
%Clear the events in which the instantaneou amplitude values are higher
%than threshold
keep_event(find(keep_event(:,1) == 0),:) = [];

%Total events
fprintf('\nNumber of events: %d\n',length(keep_event));

clear filtered_data_35_65

%% Get unique events

% for e = 1:size(pfc_data,1)
%     same_block = find(event(:,1) == e);
%
%     dif_vec = diff(event(same_block,2));
%
%     for j = 1:length(dif_vec)
%         if dif_vec(j) < 0.1*fs
%             event(same_block(j+1),:) = NaN;
%         end
%
%     end
%
% end
% event(isnan(event(:,1)),:) = [];


%% Plot events
if isempty(keep_event) == 0
    
    %     qst = input('Show the detected events? 0 - No; 1 - Yes: ');
    qst = 0;
    
    if qst == 1
        
        %Excluded events matrix
        matrix_exc = [];
        
        
        %Plot window in seconds
        window = 5;
        t = linspace(0,window*fs,window*fs);
        %         blocks_pfc = floor(matrix_pfc_filt/(window*fs)) + 1;
        
        %Concatenate the block of n seconds and the timestamps
        %         b_pfc_plot = [blocks_pfc' matrix_pfc_filt'];
        b_plot = keep_event;
        
        %Plotting loop
        for e = 1:size(b_plot,1)
            
            %Ploted block
            range = b_plot(e,1);
            
            %Range 10s
            bl = ceil(b_plot(e,2)/(window*fs));
            a = (bl-1)*(window*fs)+1;
            b = a + window*fs - 1;
            range2 = linspace(a,b,window*fs);
            
            
            %             if b_plot(e,2) >= 5*fs
            %                 if b_plot(e,4) <= 25*fs
            %                     range2 = b_plot(e,2)-5*fs+1 : b_plot(e,4) + 5*fs;
            %                 else
            %                     range2 = (block_length*fs) - 10*fs+1 : block_length * fs;
            %                 end
            %             else
            %                 range2 = 1: 10 *fs;
            %             end
            
            add_pl1 = max(pfc_data(range,range2));
            t = 1:length(range2);
            
            %Plot the period
            %PFC
            
            plot(t,pfc_data(range,range2),'Color','k')
            hold on
            plot(t,filtered_data_10_16(range,range2) + add_pl1,'Color','b')
            hold on
            plot(t,hil_filt_10_16(range,range2) + add_pl1,'Color','r')
            xlabel('Seconds')
            title('PFC')
            
            %Plot the event
            dif = (b_plot(e,4) - b_plot(e,2));
            bl = ceil(b_plot(e,2)/(window*fs));
            x1 = b_plot(e,2) - (bl-1)*(window*fs);
            x2 = dif;
            y1 = min(pfc_data(range,range2));
            y2 = max(hil_filt_10_16(range,range2)+add_pl1)-min(pfc_data(range,range2));
            position = [x1 y1 x2 y2];
            
            %Plot the detected event
            rectangle('Position',position,'EdgeColor','r')
            
            legend('PFC LFP','Filtered Data','Filtered Data Envelope','Event','location','southoutside')
            hold off
            
            %Manual correction
            [~,~,button] = ginput(1);
            
            switch button
                case 48
                    matrix_exc = [matrix_exc e];
                    %                 matrix_hipo_filt(e) = [];
            end
        end
        
        %Excluded the manually selected event
        keep_event(matrix_exc,:) = [];
        
    end
else
    disp('Events not found!')
end

timestamps = keep_event;

% sav = sprintf('/home/cleiton/Documents/ICs/Ikaro/Dados/Spindle/R%d_%d_TESTE.mat',rat,ratnum);
% save(sav,'timestamps','thresh1','thresh2')


end