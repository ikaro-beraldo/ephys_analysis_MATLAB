%  140 80 140 200 300
%   1          1.5

clear all
name = 'B2_D1';

load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data',name,'blocked_data.mat'),'LFP1','LFP3','LFP2')
load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data',name,'GMM_Classification.mat'),'GMM')
load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\ripple',name),'ripple_accepted','ripple_result_block_2')
load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\delta',name),'delta_blocks')
load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\spindle',name),'spindle_blocks')

av_ = mean(LFP3');
av_abs = mean(abs(LFP3'));
av_t = mean(LFP3','all');
av_abs_t = mean(abs(LFP3'),'all');

% Exclude av_abs > 500 from LFP
LFP3(av_abs > 500,:) = 0;

subplot(2,1,1)
plot(av_)
subplot(2,1,2)
plot(av_abs)

Check_Ripple_Delta_Spindle(GMM.All_Sort, LFP1, LFP3, ripple_accepted, delta_blocks, spindle_blocks)
%% 
% for ii = 1:size(delta_blocks,1)
%     idx = delta_blocks(ii,1);
%     a = delta_blocks(ii,2);
%     b = delta_blocks(ii,4);
% 
%     delta_amp(ii) = max(LFP1(idx,a:b));
%     delta_int(ii) = trapz(LFP1(idx,a:b));
% 
%     plot(1:10000,LFP1(idx,:))
%     hold on
%     plot(a:b,LFP1(idx,a:b))
%     text(b,delta_amp(ii),string(delta_amp(ii)))
%     hold off
%     ylim([-1500 1500])    
%     pause(0.5)
% end
% scatter(delta_int,delta_amp)

% Detect ripple

CA1 = reshape(LFP3(GMM.All_Sort==2,:)',[numel(LFP3(GMM.All_Sort==2,:)) 1])';
NREM = find(GMM.All_Sort==2);
fs = 1000;
b_length = 10;
[result] = detectRipplesOfflineKIB(CA1, fs, false, [], 4);

ripple_block = result{1};
idx = ceil((ripple_block(:,1)/fs)/b_length);
timestamps = ripple_block(:,1:3)-(idx-1)*(fs*b_length);
timestamps(timestamps(:,1) > fs*b_length,:) = [];
timestamps(timestamps(:,2) > fs*b_length,2) = 9999;
timestamps(timestamps(:,3) > fs*b_length,3) = 10000;

ripple_result_block_3 = [NREM(idx) timestamps ripple_block(:,4)];

av_ = mean(LFP3);
av_abs = mean(abs(LFP3));
av_abs_t_ratio= mean(abs(LFP3'),'all')/133;

Check_Ripple_Delta_Spindle(GMM.All_Sort, LFP1, LFP3./av_abs_t_ratio, ripple_accepted, delta_blocks, spindle_blocks)


%% Clustering ripples in 3 different clusters based on RMS and Sharp-wave amplitude

filtered_LFP2 = eegfilt2(LFP2,fs,140,220);

% Loop for each ripple
ripple_rms = nan(size(ripple_result_block_3,1),1);
ripple_amp = nan(size(ripple_result_block_3,1),1);
ripple_duration = nan(size(ripple_result_block_3,1),1);

for ii = 1:size(ripple_result_block_3,1)
    idx = ripple_result_block_3(ii,1);
    a = ripple_result_block_3(ii,2);
    b = ripple_result_block_3(ii,4);
    ripple_rms(ii) = rms(filtered_LFP2(idx,a:b)./av_abs_t_ratio);
    ripple_amp(ii) = min(LFP2(idx,a:b)./av_abs_t_ratio);
    ripple_duration(ii) = (b-a)/fs;
end

% bagulho = true;
% while bagulho
%     k_idx = kmeans([ripple_amp ripple_rms],3);
%
%     subplot(2,1,1)
%     scatter(ripple_amp,ripple_rms)
%
%     subplot(2,1,2)
%     scatter(ripple_amp(k_idx==1),ripple_rms(k_idx==1))
%     hold on
%     scatter(ripple_amp(k_idx==2),ripple_rms(k_idx==2))
%     scatter(ripple_amp(k_idx==3),ripple_rms(k_idx==3))
%     hold off
%
%     pause(1)
% end


%% Normalize

ripple_amp = zscore(ripple_amp);
ripple_rms = zscore(ripple_rms);
ripple_duration = zscore(ripple_duration);

%% GMM

bagulho = true;
while bagulho
    GMModel = fitgmdist([ripple_amp ripple_rms],3);
    pp = posterior(GMModel,[ripple_amp ripple_rms]);

    [~,gmmIDX] = max(pp,[],2);

    subplot(2,1,1)
    scatter(ripple_amp,ripple_rms)

    subplot(2,1,2)
    scatter(ripple_amp(gmmIDX==1),ripple_rms(gmmIDX==1))
    hold on
    scatter(ripple_amp(gmmIDX==2),ripple_rms(gmmIDX==2))
    scatter(ripple_amp(gmmIDX==3),ripple_rms(gmmIDX==3))
    hold off

    pause(0.2)

    % Check the user choices
    [~,~,button] = ginput(1);
    switch button
        case 29 % Left --> continues
        case 32 % Accept the upper cluster as a ripple cluster
            bagulho = false;
            [~, upper] = max([mean(ripple_rms(gmmIDX==1)) mean(ripple_rms(gmmIDX==2)) mean(ripple_rms(gmmIDX==3))]);

            ripple_sharp_wave = ripple_result_block_3(gmmIDX==upper,:);

            % Plot the new ripples
            figure
            scatter(ripple_amp(gmmIDX==upper),ripple_rms(gmmIDX==upper))
            pause
            save(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\ripple',name),"ripple_revisited",'-append')
        case 50 % Accept the upper two clusters


            bagulho = false;
            [~, upper] = sort([mean(ripple_rms(gmmIDX==1)) mean(ripple_rms(gmmIDX==2)) mean(ripple_rms(gmmIDX==3))],2,'descend');

            ripple_sharp_wave = ripple_result_block_3(gmmIDX==upper(1),:);
            ripple_no_sw = ripple_result_block_3(gmmIDX==upper(2),:);
            ripple_revisited = ripple_result_block_3(gmmIDX==upper(1) | gmmIDX==upper(2),:);

            % Plot the new ripples
            figure
            subplot(2,1,1)
            scatter(ripple_amp,ripple_rms)
            subplot(2,1,2)
            scatter(ripple_amp(gmmIDX==upper(1)),ripple_rms(gmmIDX==upper(1)))
            hold on
            scatter(ripple_amp(gmmIDX==upper(2)),ripple_rms(gmmIDX==upper(2)))
            pause
            save(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\ripple',name),"ripple_revisited",'ripple_sharp_wave','ripple_no_sw','-append')
    end
end


%% Ripple RMS/ Epoch Differentiate High Gamma Epochs

clear all
for ii = [6]
    for jj = [1]
        name = sprintf('B%d_D%d',ii,jj);
        filename = fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\ripple',name);
        load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data',name,'blocked_data.mat'),'LFP3','LFP2')        
        load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\ripple',name),'ripple_revisited')

        % Unique case
        if ii == 4
            LFP3 = LFP2;
        end

        % Filt data
        filtered_LFP3 = eegfilt2(LFP3,1000,140,220);

        % RMS total
        rms_total = rms(filtered_LFP3(ripple_revisited(:,1),:)');
        ripple_revisited_rms = nan(size(ripple_revisited,1),1);
        ripple_amp = nan(size(ripple_revisited,1),1);

        for zz = 1:size(ripple_revisited,1)

            idx = ripple_revisited(zz,1);
            a = ripple_revisited(zz,2);
            b = ripple_revisited(zz,4);

            ripple_revisited_rms(zz) = rms(filtered_LFP3(idx,a:b));
            ripple_amp(zz) = min(filtered_LFP3(idx,a:b));

        end

        histogram(ripple_revisited_rms./rms_total',100)

        bagulho = true;
        while bagulho
            % Use gmm to separate
            GMModel = fitgmdist(ripple_revisited_rms./rms_total',3);
            pp = posterior(GMModel,ripple_revisited_rms./rms_total');

            [~,gmmIDX] = max(pp,[],2);

            f = figure;
            % Plot RMS vs Amp

            subplot(3,1,2)
            scatter(ripple_amp(gmmIDX==1),ripple_revisited_rms(gmmIDX==1))
            hold on
            scatter(ripple_amp(gmmIDX==2),ripple_revisited_rms(gmmIDX==2))
            scatter(ripple_amp(gmmIDX==3),ripple_revisited_rms(gmmIDX==3))
            hold off
            xlabel('Amp')
            ylabel('RMS')

         
            % Plot RMS/RMStotal vs RMS
            ripple_rmsVSrms_total = ripple_revisited_rms./rms_total';
            subplot(3,1,3)
            scatter(ripple_rmsVSrms_total(gmmIDX==1),ripple_revisited_rms(gmmIDX==1))
            hold on
            scatter(ripple_rmsVSrms_total(gmmIDX==2),ripple_revisited_rms(gmmIDX==2))
            scatter(ripple_rmsVSrms_total(gmmIDX==3),ripple_revisited_rms(gmmIDX==3))
            hold off
            xlabel('RMS Ratio')
            ylabel('RMS')

            subplot(3,1,1)
            histogram(ripple_rmsVSrms_total(gmmIDX==1),100)
            hold on
            histogram(ripple_rmsVSrms_total(gmmIDX==2),100)
            histogram(ripple_rmsVSrms_total(gmmIDX==3),100)
            hold off
            xlabel('RMS Ratio')

            pause(0.1)

            % Check the user choices
            [~,~,button] = ginput(1);
            switch button
                case 29 % Left --> 
                close(f)

                case 48 % 0 - Do not exclude ripples
                    bagulho = false;
                    ripple_accepted = ripple_revisited;
                    ripple_excluded = [];

                case 49 % 1 - Exclude the events lower than 1.5 ratio
                    bagulho = false;
                    ripple_accepted = ripple_revisited(ripple_rmsVSrms_total >= 1,:);
                    ripple_excluded = ripple_revisited(ripple_rmsVSrms_total < 1,:);

                case 50 % 2 - Accept the upper two clusters
                    bagulho = false;
                    [~, upper] = sort([mean(ripple_rmsVSrms_total(gmmIDX==1)) mean(ripple_rmsVSrms_total(gmmIDX==2)) mean(ripple_rmsVSrms_total(gmmIDX==3))],2,'descend');

                    ripple_accepted = ripple_revisited(gmmIDX==upper(1) | gmmIDX==upper(2),:);
                    ripple_excluded = ripple_revisited(gmmIDX==upper(3),:);

            end

        end

        % Plot excluded ripples
                    figure
                    for zz = 1:size(ripple_excluded)
                        idx = ripple_excluded(zz,1);
                        a = ripple_excluded(zz,2);
                        b = ripple_excluded(zz,4);

                        subplot(2,1,1)
                        plot(1:10000,LFP3(idx,:))
                        hold on
                        % ripple
                        plot(a:b,LFP3(idx,a:b))
                        ylim([-1000 1000])
                        hold off

                        % Filtered signal

                        subplot(2,1,2)
                        plot(1:10000,filtered_LFP3(idx,:))
                        hold on
                        % ripple
                        plot(a:b,filtered_LFP3(idx,a:b))
                        ylim([-200 200])
                        hold off
                        pause(0.3)
                    end
                    pause
        % Save accepted ripples
        save(filename,'ripple_accepted','ripple_excluded','-append')
        close all
    end
end
