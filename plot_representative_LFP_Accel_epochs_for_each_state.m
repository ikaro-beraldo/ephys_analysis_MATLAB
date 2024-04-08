% Plot representative LFP
clear all

lfp_data = 'E:\Barnes Maze - Mestrad\dados matlab\blocked_data\B3_D1\blocked_data.mat';
states_data = 'E:\Barnes Maze - Mestrad\dados matlab\blocked_data\B3_D1\GMM_Classification.mat';
psd_data = 'E:\Barnes Maze - Mestrad\dados matlab\blocked_data\PSD\B3_D1.mat';
load(states_data,'GMM')
load(lfp_data,'LFP1','LFP2','Accel')
load(psd_data,'PSD')


% Info regarding the number of epochs ploted
sel_epoch = 1880;
number_epochs = 75;

f2 = figure('PaperSize',[21 11.25]);

%% Plot parameters

figure_parameters.color.lfp(4,:) = [0 0 0]; % Unclassified epochs
figure_parameters.color.lfp(3,:) = [0.9290, 0.6940, 0.1250];
figure_parameters.color.lfp(2,:)=[0 0.4470 0.7410];
figure_parameters.color.lfp(1,:) = [0.3 0.3 0.3];
figure_parameters.color.awake = [0.9290, 0.6940, 0.1250];
figure_parameters.color.nrem =[0 0.4470 0.7410];
figure_parameters.color.rem = [0.3 0.3 0.3];
figure_parameters.red = [1 0.2 0.3];

figure_parameters.fontsize=20;
figure_parameters.lw=2;

% Sizes
figure_parameters.fontsize=20;
figure_parameters.scatter_size=5;
figure_parameters.scatter_size2 = 150;
figure_parameters.edges=-3:0.1:6;
figure_parameters.lw=1.5;

% General settings
figure_parameters.transparecy_fa=.9;
figure_parameters.limx=[-1.5 8];
figure_parameters.limy=[-3 10];

%% Spectrogram specs

FS = 1000;
block_length = 10; % Seconds

%PARAMS
WindowLength = 0.25; %em  segundos
WindowLength = WindowLength*FS;
Overlap = 0.5*FS;
NFFT = 2^15;

%% Define the best epochs for each state

selected_epoch_vector = sel_epoch-number_epochs:sel_epoch+number_epochs;

NREM = find(GMM.All_Sort == 2);
REM = find(GMM.All_Sort == 1);
WK = find(GMM.All_Sort == 3);

[~,NREM_s] = sort(GMM.LFP_used(intersect(NREM,selected_epoch_vector)));
[~,REM_s] = sort(GMM.LFP_used(intersect(REM,selected_epoch_vector)));
[~,WK_s] = sort(GMM.EMG_used(intersect(WK,selected_epoch_vector)));


%% 1 - mPFC LFP plot

time_vector = linspace(0,block_length,FS*block_length);

counter = 1;
for ii = [REM(REM_s(end-1)) NREM(NREM_s(1)) WK(WK_s(end-1))]
    figure
    % Plot the LFP
    subplot(3,1,1)
    plot(time_vector,LFP1(ii,:),'Color',figure_parameters.color.lfp(counter,:))
    ylim([-1000 1000])
    xlim([0 time_vector(end)])
    set(gca,'fontsize',figure_parameters.fontsize)
    set(gca,'Linewidth',figure_parameters.lw)
    set(gca,'Tickdir','out')
    ylabel('mV')
    xticklabels({})
    box off
    title(sprintf('%d',ii))

    subplot(3,1,2)
    plot(time_vector,LFP2(ii,:),'Color',figure_parameters.color.lfp(counter,:))
    ylim([-1000 1000])
    xlim([0 time_vector(end)])
    set(gca,'fontsize',figure_parameters.fontsize)
    set(gca,'Linewidth',figure_parameters.lw)
    set(gca,'Tickdir','out')
    ylabel('mV')
    xticklabels({})
    box off

    subplot(3,1,3)
    plot(time_vector,Accel(ii,:),'Color',figure_parameters.color.lfp(counter,:))
    ylim([-0.5 0.5])
    xlim([0 time_vector(end)])
    set(gca,'fontsize',figure_parameters.fontsize)
    set(gca,'Linewidth',figure_parameters.lw)
    set(gca,'Tickdir','out')
    ylabel('mV')
    xticklabels({})
    box off

    counter = counter + 1;

end

%% Plot representative PSDs

counter = 1;
for ii = [REM(REM_s(end-1)) NREM(NREM_s(1)) WK(WK_s(end-1))]
    figure
    % Plot the LFP
    subplot(1,2,1)
    loglog(PSD.auxF,PSD.LFP1(:,ii),'Color',figure_parameters.color.lfp(counter,:))
    xlim([0 100])
    set(gca,'fontsize',figure_parameters.fontsize)
    set(gca,'Linewidth',figure_parameters.lw)
    set(gca,'Tickdir','out')
    ylabel('Power')
    xlabel('Frequency')
    box off
    title(sprintf('%d',ii))

    subplot(1,2,2)
    loglog(PSD.auxF,PSD.LFP2(:,ii),'Color',figure_parameters.color.lfp(counter,:))
    xlim([0 100])
    set(gca,'fontsize',figure_parameters.fontsize)
    set(gca,'Linewidth',figure_parameters.lw)
    set(gca,'Tickdir','out')
    ylabel('Power')
    xlabel('Frequency')
    box off

%     subplot(3,1,3)
%     plot(time_vector,Accel(ii,:),'Color',figure_parameters.color.lfp(counter,:))
%     ylim([-0.5 0.5])
%     xlim([0 time_vector(end)])
%     set(gca,'fontsize',figure_parameters.fontsize)
%     set(gca,'Linewidth',figure_parameters.lw)
%     set(gca,'Tickdir','out')
%     ylabel('mV')
%     xticklabels({})
%     box off

    counter = counter + 1;

end

%% Plot all states (mpfc and Ca1 on the same plot)

REM_p = REM(REM_s(end-1));
NREM_p = NREM(NREM_s(1));
WK_p = WK(WK_s(end-1));

figure
% Plot REM
subplot(1,3,1)
loglog(PSD.auxF,PSD.LFP1(:,REM_p),'--','Color',figure_parameters.color.lfp(1,:),'LineWidth',figure_parameters.lw)
hold on
loglog(PSD.auxF,PSD.LFP2(:,REM_p),'Color',figure_parameters.color.lfp(1,:),'LineWidth',figure_parameters.lw)
xlim([0 100])
set(gca,'fontsize',figure_parameters.fontsize)
set(gca,'Linewidth',figure_parameters.lw)
set(gca,'Tickdir','out')
ylabel('Power')
xlabel('Frequency (Hz)')
box off
title(sprintf('%d',REM_p))

% Plot NREM
subplot(1,3,2)
loglog(PSD.auxF,PSD.LFP1(:,NREM_p),'--','Color',figure_parameters.color.lfp(2,:),'LineWidth',figure_parameters.lw)
hold on
loglog(PSD.auxF,PSD.LFP2(:,NREM_p),'Color',figure_parameters.color.lfp(2,:),'LineWidth',figure_parameters.lw)
xlim([0 100])
set(gca,'fontsize',figure_parameters.fontsize)
set(gca,'Linewidth',figure_parameters.lw)
set(gca,'Tickdir','out')
ylabel('Power')
xlabel('Frequency (Hz)')
box off
title(sprintf('%d',NREM_p))

% Plot WK
subplot(1,3,3)
loglog(PSD.auxF,PSD.LFP1(:,WK_p),'--','Color',figure_parameters.color.lfp(3,:),'LineWidth',figure_parameters.lw)
hold on
loglog(PSD.auxF,PSD.LFP2(:,WK_p),'Color',figure_parameters.color.lfp(3,:),'LineWidth',figure_parameters.lw)
xlim([0 100])
set(gca,'fontsize',figure_parameters.fontsize)
set(gca,'Linewidth',figure_parameters.lw)
set(gca,'Tickdir','out')
ylabel('Power')
xlabel('Frequency (Hz)')
box off
title(sprintf('%d',WK_p))


    %% Save it
% 
% set(gcf,'color','white')
% set(f2,'PaperPositionMode','auto')
% % set(f2,'PaperOrientation','landscape')
% sgtitle(sprintf('Representative epochs around %d',sel_epoch),'fontsize',figure_parameters.fontsize*2.2)
% print('-fillpage',sprintf('Representative epochs around %d',sel_epoch),'-dpdf','-r0',f2,'-vector')
% 
% pause(10)