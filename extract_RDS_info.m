function extract_RDS_info(directory,names)

for ii = 1:length(names)
    disp(names{ii})
    clearvars -except directory names counter ii s_p
    counter = 1;

    %% RIPPLE

    % Complete filename
    filename1 = fullfile(directory,names{ii},'blocked_data.mat');
    filename2 = fullfile(directory,names{ii},'GMM_Classification.mat');
    filename3 = fullfile(directory,'ripple',names{ii});
    load(filename1,'LFP2','LFP3','fs')
    load(filename2)
    load(filename3,'ripple_result_block_2','ripple_result_block_3','ripple_accepted')

    % Calculate duration
    duration.ripple_2 = (ripple_result_block_2(:,4) - ripple_result_block_2(:,2))/fs;
    duration.ripple_3 = (ripple_result_block_3(:,4) - ripple_result_block_3(:,2))/fs;
    duration.ripple = (ripple_accepted(:,4) - ripple_accepted(:,2))/fs;

   
    % Get the peak frequency LFP2
    l = 140;
    h = 220;
    [filtered_data]=eegfilt2(LFP2,fs,l,[]);
    [filtered_data]=eegfilt2(filtered_data,fs,[],h);

     % Calculate the RMS
    for rr = 1:size(ripple_result_block_2,1)
        a = ripple_result_block_2(rr,2);
        b = ripple_result_block_2(rr,4);
        idx = ripple_result_block_2(rr,1);
        RMS.ripple_2(rr) = rms(filtered_data(idx,a:b));
        % Calculate PEAK
        PEAK.ripple_2(rr) = max(abs(filtered_data(idx,a:b)));    
    end
    
    %fft_test
    L = 10000;
    f = fs*(0:(L/2))/L;
    % Loop to compute FFT
    for rr = 1:size(ripple_result_block_2,1)
        a = ripple_result_block_2(rr,2);
        b = ripple_result_block_2(rr,4);
        idx = ripple_result_block_2(rr,1);
        Y = fft(filtered_data(idx,a:b),10000);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        [~,m] = max(P1);    % Get idx from the frequency with highest power
        FREQ.ripple_2(rr) = f(m); % Get the ripple frequency
    end

    % Get the peak frequency LFP3
    l = 140;
    h = 220;
    [filtered_data]=eegfilt2(LFP3,fs,l,[]);
    [filtered_data]=eegfilt2(filtered_data,fs,[],h);

    for rr = 1:size(ripple_result_block_3,1)
        a = ripple_result_block_3(rr,2);
        b = ripple_result_block_3(rr,4);
        idx = ripple_result_block_3(rr,1);
        RMS.ripple_3(rr) = rms(filtered_data(idx,a:b));
        % Calculate PEAK
        PEAK.ripple_3(rr) = max(abs(filtered_data(idx,a:b)));
    end

    %fft_test
    L = 10000;
    f = fs*(0:(L/2))/L;
    % Loop to compute FFT
    for rr = 1:size(ripple_result_block_3,1)
        a = ripple_result_block_3(rr,2);
        b = ripple_result_block_3(rr,4);
        idx = ripple_result_block_3(rr,1);
        Y = fft(filtered_data(idx,a:b),10000);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        [~,m] = max(P1);    % Get idx from the frequency with highest power
        FREQ.ripple_3(rr) = f(m); % Get the ripple frequency
    end
    
    % Get accepted ripple RMS
    l = 140;
    h = 220;
    if strcmp(names{ii}(1:2),'B4')
        % Get the peak frequency LFP3
        [filtered_data]=eegfilt2(LFP2,fs,l,[]);
    else
        % Get the peak frequency LFP3
        [filtered_data]=eegfilt2(LFP3,fs,l,[]);
    end
    [filtered_data]=eegfilt2(filtered_data,fs,[],h);
    for rr = 1:size(ripple_accepted,1)
        a = ripple_accepted(rr,2);
        b = ripple_accepted(rr,4);
        idx = ripple_accepted(rr,1);
        RMS.ripple(rr) = rms(filtered_data(idx,a:b));        
        % Calculate PEAK
        PEAK.ripple(rr) = max(abs(filtered_data(idx,a:b)));
    end

    %FFT_test
    L = 10000;
    f = fs*(0:(L/2))/L;
    % Loop to compute FFT
    for rr = 1:size(ripple_accepted,1)
        a = ripple_accepted(rr,2);
        b = ripple_accepted(rr,4);
        idx = ripple_accepted(rr,1);
        Y = fft(filtered_data(idx,a:b),10000);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        [~,m] = max(P1);    % Get idx from the frequency with highest power
        FREQ.ripple(rr) = f(m); % Get the ripple frequency
    end

    sav = fullfile(directory,'ripple',names{ii});
    save(sav,'duration','RMS','PEAK','FREQ','-append')
    clear duration RMS PEAK FREQ

    %% SPINDLE

    % Complete filename
    filename1 = fullfile(directory,names{ii},'blocked_data.mat');
    filename2 = fullfile(directory,names{ii},'GMM_Classification.mat');
    filename3 = fullfile(directory,'spindle',names{ii});
    load(filename1,'LFP1','fs')
    load(filename2)
    load(filename3,'spindle_blocks')

    % Calculate duration
    duration.spindle = (spindle_blocks(:,4) - spindle_blocks(:,2))/fs;
    
    % Calculate the RMS
    for rr = 1:size(spindle_blocks,1)
        a = spindle_blocks(rr,2);
        b = spindle_blocks(rr,4);
        idx = spindle_blocks(rr,1);
        RMS.spindle(rr) = rms(LFP1(idx,a:b));
    end

    % Calculate PEAK
    for rr = 1:size(spindle_blocks,1)
        PEAK.spindle(rr) = abs(LFP1(spindle_blocks(rr,1),spindle_blocks(rr,3)));
    end

     % Get the peak frequency LFP3
    l = 10;
    h = 16;
    [filtered_data]=eegfilt2(LFP1,fs,l,[]);
    [filtered_data]=eegfilt2(filtered_data,fs,[],h);

    %fft_test
    L = 10000;
    f = fs*(0:(L/2))/L;
    % Loop to compute FFT
    for rr = 1:size(spindle_blocks,1)
        a = spindle_blocks(rr,2);
        b = spindle_blocks(rr,4);
        idx = spindle_blocks(rr,1);
        Y = fft(filtered_data(idx,a:b),10000);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        [~,m] = max(P1);    % Get idx from the frequency with highest power
        FREQ.spindle(rr) = f(m); % Get the ripple frequency
    end
    
    sav = fullfile(directory,'spindle',names{ii});
    save(sav,'duration','RMS','PEAK','FREQ','-append')

end
end
