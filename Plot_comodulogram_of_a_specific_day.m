% Plot representative data 
clear all

ca1_comodulation = 'E:\Barnes Maze - Mestrad\dados matlab\comodulation\CA1\B3_D1.mat';
ca1_mpfc_comodulation = 'E:\Barnes Maze - Mestrad\dados matlab\comodulation\CA1-mPFC\B3_D1.mat';
states_data = 'E:\Barnes Maze - Mestrad\dados matlab\blocked_data\B3_D1\GMM_Classification.mat';
load(states_data,'GMM')

% Get the states
NREM = find(GMM.All_Sort == 2);
REM = find(GMM.All_Sort == 1);
WK = find(GMM.All_Sort == 3);

% Info regarding the number of epochs ploted
sel_epoch = 1880;
number_epochs = 75;
time_vector = sel_epoch - number_epochs: sel_epoch + number_epochs;

NREM_s = intersect(NREM,time_vector);
REM_s = intersect(REM,time_vector);
WK_s = intersect(WK,time_vector);

NREM_s = NREM;
REM_s = REM;
WK_s = WK;

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


% CA1 - 1 comodulation    
sp = subplot(1,3,1);
load(ca1_comodulation,'Comodulogram_all','parameters')
matrix_plot = mean(Comodulogram_all(:,:,REM_s),3)';
contourf(parameters.PhaseFreqVector+parameters.PhaseFreq_BandWidth/2, parameters.AmpFreqVector+parameters.AmpFreq_BandWidth/2, matrix_plot,120,'lines','none')
ylabel('Amp. Frequency (Hz)' )
% xticks(position_idx)
set(gca,'fontsize',figure_parameters.fontsize)
set(gca,'Linewidth',figure_parameters.lw)
set(gca,'Tickdir','out')
sp.Position
colorbar
clim([0.0006 0.0018])
axis xy
colormap(magma)
t=xlabel('Phase Frequency (Hz)');
t.Color = 'k';
sp.Position

sp = subplot(1,3,2);
% NREM
matrix_plot = mean(Comodulogram_all(:,:,NREM_s),3)';
contourf(parameters.PhaseFreqVector+parameters.PhaseFreq_BandWidth/2, parameters.AmpFreqVector+parameters.AmpFreq_BandWidth/2, matrix_plot,120,'lines','none')
ylabel('Amp. Frequency (Hz)' )
% xticks(position_idx)
set(gca,'fontsize',figure_parameters.fontsize)
set(gca,'Linewidth',figure_parameters.lw)
set(gca,'Tickdir','out')
sp.Position
colorbar
clim([0.001 0.0025])
axis xy
colormap(magma)
t=xlabel('Phase Frequency (Hz)');
t.Color = 'k';
sp.Position

%WK
sp = subplot(1,3,3);
matrix_plot = mean(Comodulogram_all(:,:,WK_s),3)';
contourf(parameters.PhaseFreqVector+parameters.PhaseFreq_BandWidth/2, parameters.AmpFreqVector+parameters.AmpFreq_BandWidth/2, matrix_plot,120,'lines','none')
% imagesc(parameters.PhaseFreqVector,parameters.AmpFreqVector,smoothn(matrix_plot,10));
ylabel('Amp. Frequency (Hz)' )
% xticks(position_idx)
set(gca,'fontsize',figure_parameters.fontsize)
set(gca,'Linewidth',figure_parameters.lw)
set(gca,'Tickdir','out')
sp.Position
colorbar
clim([0.001 0.0018])
axis xy
colormap(magma)
t=xlabel('Phase Frequency (Hz)');
t.Color = 'k';
sp.Position


%% Print PDF

f2.Renderer='Painters';
set(gcf,'color','white')
set(f2,'PaperPositionMode','auto')
% set(f2,'PaperOrientation','landscape')
sgtitle(sprintf('Comodulogram each state representative %d',sel_epoch),'fontsize',figure_parameters.fontsize*2.2)
print('-fillpage',sprintf('Comodulogram each state representative 2 %d',sel_epoch),'-dpdf','-r0',f2,'-vector')

% pause(10)