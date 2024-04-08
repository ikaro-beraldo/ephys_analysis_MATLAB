% Plot representative data 
clear all

spectral_data = 'E:\Barnes Maze - Mestrad\dados matlab\blocked_data\spectral_coherence13\B3_D1.mat';
phase_data = 'E:\Barnes Maze - Mestrad\dados matlab\blocked_data\phase_coherence13\B3_D1.mat';
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


% Phase Coherence    
sp = subplot(3,1,1);
time_vector = sel_epoch - number_epochs: sel_epoch + number_epochs;
time_vector_min = time_vector*block_length/60;
load(phase_data,'Call13','f')
imagesc(time_vector_min,f,smoothn(Call13(:,time_vector),10));

ylabel('Frequency (Hz)' )
xlim([time_vector_min(1) time_vector_min(end)])
% xticks(position_idx)
set(gca,'fontsize',figure_parameters.fontsize)
set(gca,'Linewidth',figure_parameters.lw)
set(gca,'Tickdir','out')
sp.Position
colorbar
clim([0.3 0.8])
axis xy
ylim([0 30])
colormap(magma)
t=xlabel('Minutes');
t.Color = 'k';
sp.Position

%% Plot Spectral coherence

sp = subplot(3,1,2);
time_vector = sel_epoch - number_epochs: sel_epoch + number_epochs;
time_vector_min = time_vector*block_length/60;
load(spectral_data,'Cxy_all13','F')
imagesc(time_vector_min,F,smoothn(Cxy_all13(:,time_vector),10));

ylabel('Frequency (Hz)' )
xlim([time_vector_min(1) time_vector_min(end)])
% xticks(position_idx)
set(gca,'fontsize',figure_parameters.fontsize)
set(gca,'Linewidth',figure_parameters.lw)
set(gca,'Tickdir','out')
sp.Position
colorbar
clim([0.2 0.7])
axis xy
ylim([0 30])
colormap(magma)
t=xlabel('Minutes');
t.Color = 'k';
sp.Position

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
sgtitle(sprintf('Epoch %d Phase and Spectral coherence',sel_epoch),'fontsize',figure_parameters.fontsize*2.2)
print('-fillpage',sprintf('Epoch %d Phase and Spectral coherence',sel_epoch),'-dpdf','-r0',f2,'-vector')

% pause(10)