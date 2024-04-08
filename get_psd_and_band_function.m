%% Function to get the PSD and Bands

function get_psd_and_band_function(directory,names)

for ii = 1:length(names)
    
    mkdir(fullfile(directory,'PSD'))
    clearvars -except directory names counter ii
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
    
    %% 1 - Pre-processing
    s_p = Signal_Processing;
    
    LFP1 = s_p.basic_preprocessing(LFP1,fs,[1 300],false);
    LFP2 = s_p.basic_preprocessing(LFP2,fs,[1 300],false);
    LFP3 = s_p.basic_preprocessing(LFP3,fs,[1 300],false);
    
    %% PSD - Pwelch
    PSD = struct;
    
    [PSD.LFP1, PSD.auxF] = s_p.calc_PSD(LFP1, fs);
    [PSD.LFP2, PSD.auxF] = s_p.calc_PSD(LFP2, fs);
    [PSD.LFP3, PSD.auxF] = s_p.calc_PSD(LFP3, fs);
    
    PSD.LFP1 = PSD.LFP1';
    PSD.LFP2 = PSD.LFP2';
    PSD.LFP3 = PSD.LFP3';
    
    save_psd = fullfile(directory,'PSD',sprintf('%s',names{ii}));
    save(save_psd,'PSD','fs')
    
    %% Get frequency bands
    
    f = PSD.auxF >= 1 & PSD.auxF <= 4;
    BANDS.Delta.LFP1 = nansum(PSD.LFP1(f,:));
    BANDS.Delta.LFP2 = nansum(PSD.LFP2(f,:));
    BANDS.Delta.LFP3 = nansum(PSD.LFP3(f,:));
    
    f = PSD.auxF >= 5 & PSD.auxF <= 12;
    BANDS.Theta.LFP1 = nansum(PSD.LFP1(f,:));
    BANDS.Theta.LFP2 = nansum(PSD.LFP2(f,:));
    BANDS.Theta.LFP3 = nansum(PSD.LFP3(f,:));
    
    f = PSD.auxF >= 4 & PSD.auxF <= 7;
    BANDS.Theta2.LFP1 = nansum(PSD.LFP1(f,:));
    BANDS.Theta2.LFP2 = nansum(PSD.LFP2(f,:));
    BANDS.Theta2.LFP3 = nansum(PSD.LFP3(f,:));
    
    f = PSD.auxF >= 8 & PSD.auxF <= 12;
    BANDS.Theta1.LFP1 = nansum(PSD.LFP1(f,:));
    BANDS.Theta1.LFP2 = nansum(PSD.LFP2(f,:));
    BANDS.Theta1.LFP3 = nansum(PSD.LFP3(f,:));
    
    f = PSD.auxF >= 30 & PSD.auxF <= 80;
    BANDS.Gamma.LFP1 = nansum(PSD.LFP1(f,:));
    BANDS.Gamma.LFP2 = nansum(PSD.LFP2(f,:));
    BANDS.Gamma.LFP3 = nansum(PSD.LFP3(f,:));
    
    save_psd = fullfile(directory,'PSD',sprintf('%s',names{ii}));
    save(save_psd,'BANDS','-append')
    clear BANDS
    
    %% Plot pwelch for all days
    
    colors = get(0, 'defaultAxesColorOrder');
    f=figure('PaperSize', [29.7 21],'visible','off');
    
    % Restrict the plot to 200 hz
    Fidx=find(PSD.auxF<200);
    
    % mPFC
    subplot(1,3,1)
    loglog(PSD.auxF(Fidx),nanmean(PSD.LFP1(Fidx,GMM_NREM_All_Sort==1),2),'Color',colors(1,:),'linewidth',3); xlim([0.5,100]);
    hold on
    loglog(PSD.auxF(Fidx),nanmean(PSD.LFP1(Fidx,GMM_REM_All_Sort==1),2),'Color',[0.3 0.3 0.3],'linewidth',3); xlim([0.5,100]);
    loglog(PSD.auxF(Fidx),nanmean(PSD.LFP1(Fidx,GMM_WK_All_Sort==1),2),'Color',colors(3,:),'linewidth',3); xlim([0.5,100]);
    
    legend('NREM','REM','WK','FontSize',14)
    legend boxoff
    title('mPFC')
    ylim([0.000001 0.1])
    
    
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Log10 (Power)','FontSize',16)
    hold off
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    % CA1 - 1
    subplot(1,3,2)
    loglog(PSD.auxF(Fidx),nanmean(PSD.LFP2(Fidx,GMM_NREM_All_Sort==1),2),'Color',colors(1,:),'linewidth',3); xlim([0.5,100]);
    hold on
    loglog(PSD.auxF(Fidx),nanmean(PSD.LFP2(Fidx,GMM_REM_All_Sort==1),2),'Color',[0.3 0.3 0.3],'linewidth',3); xlim([0.5,100]);
    loglog(PSD.auxF(Fidx),nanmean(PSD.LFP2(Fidx,GMM_WK_All_Sort==1),2),'Color',colors(3,:),'linewidth',3); xlim([0.5,100]);
    
    legend('NREM','REM','WK','FontSize',14)
    legend boxoff
    ylim([0.000001 0.1])
    
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Log10 (Power)','FontSize',16)
    title('CA1 1')
    hold off
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    % CA1 - 2
    subplot(1,3,3)
    loglog(PSD.auxF(Fidx),nanmean(PSD.LFP3(Fidx,GMM_NREM_All_Sort==1),2),'Color',colors(1,:),'linewidth',3); xlim([0.5,100]);
    hold on
    loglog(PSD.auxF(Fidx),nanmean(PSD.LFP3(Fidx,GMM_REM_All_Sort==1),2),'Color',[0.3 0.3 0.3],'linewidth',3); xlim([0.5,100]);
    loglog(PSD.auxF(Fidx),nanmean(PSD.LFP3(Fidx,GMM_WK_All_Sort==1),2),'Color',colors(3,:),'linewidth',3); xlim([0.5,100]);
    
    legend('NREM','REM','WK','FontSize',14)
    legend boxoff
    title('CA1 2')
    ylim([0.000001 0.1])
    
    
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Log10 (Power)','FontSize',16)
    hold off
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    sgtitle(sprintf('PSD %s',names{ii}))
    print('-fillpage',fullfile(directory,'PSD',sprintf('PSD %s',names{ii})),'-dpdf','-r0',f)
    %saveas(f,fullfile(directory,sprintf('PSD %s.png',names{ii})))
    %set(gcf,'position',[300 400 1400 800]);
    
end
end
