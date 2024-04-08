% Plot representative data 
clear all

lfp_data = 'E:\Barnes Maze - Mestrad\dados matlab\blocked_data\B3_D1\blocked_data.mat';
states_data = 'E:\Barnes Maze - Mestrad\dados matlab\blocked_data\B3_D1\GMM_Classification.mat';
load(states_data,'GMM')


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
figure_parameters.lw=2;

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

%% 1 - mPFC LFP plot                                                

% 
% % Plot the LFP
% plot(time_vector,mPFC_DATA,'Color',figure_parameters.color.lfp(classification_mPFC(sp),:))
% 
% 
% hold off
% ylim([-1 1])
% xlim([0 time_vector(end,end)])
% set(gca,'fontsize',figure_parameters.fontsize)
% set(gca,'Linewidth',figure_parameters.lw)
% set(gca,'Tickdir','out')
% ylabel('mV')
% xticklabels({})

load(lfp_data,'LFP1')
LFP1(GMM.All_Sort == -1,:) = 0;
% mPFC SPECTOGRAM                                           
sp = subplot(3,1,1);

[P, FA, TA, teste] = spectrogram(reshape(LFP1(sel_epoch-number_epochs:sel_epoch+number_epochs,:)',FS*(number_epochs*2+1)*block_length,[])', WindowLength,[],linspace(0,20,1500), FS); %mPFC
clear LFP1
TA = TA/60;
imagesc(TA,FA,smoothn(abs(P),100));

ylabel('Frequency (Hz)' )
xlim([TA(1) TA(end)])
% xticks(position_idx)
set(gca,'fontsize',figure_parameters.fontsize)
set(gca,'Linewidth',figure_parameters.lw)
set(gca,'Tickdir','out')
sp.Position
colorbar
clim([3000,12000])
axis xy
ylim([0.5 20])
colormap(magma)
t=xlabel('10 seconds epochs');
t.Color = 'k';
sp.Position


clear TA PA FA


%% Plot CA1


% % CA1 LFP plot                                                PLOT 3
% subplot(4,length(position_idx),[1:length(position_idx)]+length(position_idx)*2)
% hold on
% for sp = 1:length(position_idx)
%     % Plot the LFP
%     plot(time_vector(sp,:),LFP2(sel_epoch+position_idx(sp),:),'Color',figure_parameters.color.lfp(classification_CA1(sp),:))
%     % Plot EMG 
%     plot(time_vector(sp,:),CA1_DATA.EMG_epochs(sel_epoch+position_idx(sp),:)-0.5,'Color',[0 0 0])
%     % Plot a vertical line separating the different epochs
%     xline(time_vector(sp,end))
%     % Check if it is the central epoch (position_idx == 0)
%     if position_idx(sp) == 0
%         % Plot a red line at the beginning
%         xline(time_vector(sp,1),'Color','red','LineWidth',2)
%         % Plot a red line after
%         xline(time_vector(sp,end),'Color','red','LineWidth',2)
%     end
% end
% hold off
% ylim([-1 1])
% xlim([0 time_vector(end,end)])
% set(gca,'fontsize',figure_parameters.fontsize)
% set(gca,'Linewidth',figure_parameters.lw)
% set(gca,'Tickdir','out')
% ylabel('mV')
% xticklabels({})

load(lfp_data,'LFP2')
LFP2(GMM.All_Sort == -1,:) = 0;

% CA1 SPECTOGRAM                                           % PLOT 4
sp = subplot(3,1,2);

[P, FA, TA, ~] = spectrogram(reshape(LFP2(sel_epoch-number_epochs:sel_epoch+number_epochs,:)',FS*(number_epochs*2+1)*block_length,[])', WindowLength,[],linspace(0,20,1500), FS); %CA1
clear LFP2
TA = TA/60;
% imagesc(TA,FA,PA);
imagesc(TA,FA,smoothn(abs(P),100));

ylabel('Frequency (Hz)' )
xlim([TA(1) TA(end)])
% xticks(position_idx)
set(gca,'fontsize',figure_parameters.fontsize)
set(gca,'Linewidth',figure_parameters.lw)
set(gca,'Tickdir','out')
sp.Position

colorbar
clim([5000,10000])
axis xy
ylim([0.5 20])
colormap(magma)
t=xlabel('10 seconds epochs');
t.Color = 'k';

sp.Position

clear TA PA FA


%% Plot Hipnogram
sp = subplot(3,1,3);
GMM.All_Sort(GMM.All_Sort == -1) = 0;
Hypnogram_Simplified(sp,GMM.All_Sort(sel_epoch-number_epochs:sel_epoch+number_epochs),block_length,figure_parameters.color.lfp)

sp.Position


% Print PDF

f2.Renderer='Painters';
set(gcf,'color','white')
set(f2,'PaperPositionMode','auto')
% set(f2,'PaperOrientation','landscape')
sgtitle(sprintf('Epoch %d Spectrogram',sel_epoch),'fontsize',figure_parameters.fontsize*2.2)
print('-fillpage',sprintf('Epoch %d Spectrogram',sel_epoch),'-dpdf','-r0',f2,'-vector')

pause(10)