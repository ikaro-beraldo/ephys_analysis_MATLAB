% Plot representive ripple delta and spindles
clear all

data = 'B3_D2';


% Spectrogram parameters
%PARAMS
WindowLength = 0.25; %em  segundos
FS = 1000;
WindowLength = WindowLength*FS;
Overlap = 0.5*FS;
NFFT = 2^15;


lfp_data = fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data',data,'blocked_data.mat');
states_data = fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data',data,'GMM_Classification.mat');
rds_data = fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\RDS',data);
load(states_data,'GMM')
load(lfp_data,'LFP1','LFP3')
load(rds_data,'linear')

LFP1 = reshape(LFP1',numel(LFP1),[])';
LFP3 = reshape(LFP3',numel(LFP1),[])';

LFP2_filt = eegfilt2(LFP3,1000,140,220);

ripple = linear.ripple_delta_spindle_timestamps.ripple;
delta = linear.ripple_delta_spindle_timestamps.delta;
spindle = linear.ripple_delta_spindle_timestamps.spindle;


t_vector = linspace(0,10,10000);

for ii = 30:size(ripple,1)
    figure(1)
    % mPFC
    subplot(2,1,1)
    t_vector = ripple(ii,1) - 1500:spindle(ii,3) + 1500;
    center = round(mean(t_vector));
    t_vector = center - 5000:center + 5000;
    plot(t_vector, LFP1(t_vector))
    hold on
    plot(delta(ii,1):delta(ii,3), LFP1(delta(ii,1):delta(ii,3)))
    hold on
    plot(spindle(ii,1):spindle(ii,3), LFP1(spindle(ii,1):spindle(ii,3)))
    hold off
    ylim([-1200 1200])
    xlim([t_vector(1) t_vector(end)])
%     set(gca,'fontsize',figure_parameters.fontsize)
%     set(gca,'Linewidth',figure_parameters.lw)
    set(gca,'Tickdir','out')
    ylabel('uV')
    box off

    subplot(2,1,2)    
    plot(t_vector, LFP3(t_vector))
    hold on
    plot(t_vector, -1000+LFP2_filt(t_vector))
    plot(ripple(ii,1):ripple(ii,3), LFP3(ripple(ii,1):ripple(ii,3)))
    ylim([-1200 1200])
    xlim([t_vector(1) t_vector(end)])
%     set(gca,'fontsize',figure_parameters.fontsize)
%     set(gca,'Linewidth',figure_parameters.lw)
    set(gca,'Tickdir','out')
    ylabel('uV')
    box off
    
    hold off

    % Spectrograms
    figure(2)
    subplot(2,1,1)

    [P, FA, TA, teste] = spectrogram(LFP1(t_vector), WindowLength,[],linspace(0,20,3000), FS); %mPFC
    TA = TA/60;
    imagesc(TA,FA,smoothn(abs(P),30));
    ylabel('Frequency (Hz)' )
    xlim([TA(1) TA(end)])
    % xticks(position_idx)
    set(gca,'Tickdir','out')
    colorbar
    clim([5000,30000])
    axis xy
    ylim([0.5 20])
    colormap(magma)
    t=xlabel('10 seconds epochs');
    t.Color = 'k';

    subplot(2,1,2)

    [P, FA, TA, teste] = spectrogram(LFP3(t_vector), WindowLength,[],linspace(130,250,3000), FS); %mPFC
    TA = TA/60;
    imagesc(TA,FA,smoothn(abs(P),30));
    ylabel('Frequency (Hz)' )
    xlim([TA(1) TA(end)])
    % xticks(position_idx)
    set(gca,'Tickdir','out')
    colorbar
    clim([300,500])
    axis xy
    ylim([130 250])
    colormap(magma)
    t=xlabel('10 seconds epochs');
    t.Color = 'k';

pause
end

