%% Plot PSD for each DAY
function plot_psds_function(directory,names)

counter_loop = 0;

for ii = 1:size(names,1)
    for jj = 1:size(names,2)
        
        % Load files
        open_PSD = fullfile(directory,'PSD',char(names{ii,jj}));
        load(open_PSD,'PSD')
        open_class = fullfile(directory,char(names{ii,jj}),'GMM_Classification.mat');
        load(open_class,'GMM','GMM_NREM_All_Sort','GMM_REM_All_Sort','GMM_WK_All_Sort')
        
        % Fix the classification
        GMM_NREM_All_Sort(:) = false;
        GMM_NREM_All_Sort(GMM.All_Sort == 2) = true;
        GMM_REM_All_Sort(:) = false;
        GMM_REM_All_Sort(GMM.All_Sort == 1) = true;
        GMM_WK_All_Sort(:) = false;
        GMM_WK_All_Sort(GMM.All_Sort == 3) = true;
        
        % Extract the averages
        pSD.nrem.LFP1(:,jj) = nanmean(PSD.LFP1(:,GMM_NREM_All_Sort==1),2);
        pSD.rem.LFP1(:,jj) = nanmean(PSD.LFP1(:,GMM_REM_All_Sort==1),2);
        pSD.wk.LFP1(:,jj) = nanmean(PSD.LFP1(:,GMM_WK_All_Sort==1),2);
        
        pSD.nrem.LFP2(:,jj) = nanmean(PSD.LFP2(:,GMM_NREM_All_Sort==1),2);
        pSD.rem.LFP2(:,jj) = nanmean(PSD.LFP2(:,GMM_REM_All_Sort==1),2);
        pSD.wk.LFP2(:,jj) = nanmean(PSD.LFP2(:,GMM_WK_All_Sort==1),2);
        
        pSD.nrem.LFP3(:,jj) = nanmean(PSD.LFP3(:,GMM_NREM_All_Sort==1),2);
        pSD.rem.LFP3(:,jj) = nanmean(PSD.LFP3(:,GMM_REM_All_Sort==1),2);
        pSD.wk.LFP3(:,jj) = nanmean(PSD.LFP3(:,GMM_WK_All_Sort==1),2);
        
        clear Cxy_all Call
    end
    
    counter_loop = counter_loop + 1;
    if counter_loop == 1
        colors = get(0, 'defaultAxesColorOrder');
        ff=figure('PaperSize', [29.7 21]);
        counter = 0;
    end
    
    %% mPFC
    % Phase NREM
    subplot(3,3,1)
    for pp = 1:length(pSD.nrem.LFP1(1,:))
        loglog(PSD.auxF,pSD.nrem.LFP1(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     loglog(PSD.auxF,pSD.nrem.LFP1(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     loglog(PSD.auxF,pSD.nrem.LFP1(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.nrem.LFP1(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.nrem.LFP1(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.nrem.LFP1(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Power.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('mPFC - NREM')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    ylim([0.000001 0.1])
    xlim([0.5,100])
    
    % pSD REM
    subplot(3,3,2)
    for pp = 1:length(pSD.rem.LFP1(1,:))
        loglog(PSD.auxF,pSD.rem.LFP1(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     loglog(PSD.auxF,pSD.rem.LFP1(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     loglog(PSD.auxF,pSD.rem.LFP1(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.rem.LFP1(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.rem.LFP1(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.rem.LFP1(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Power.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('mPFC - REM')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    ylim([0.000001 0.1])
    xlim([0.5,100])
    
    % Phase WK
    subplot(3,3,3)
    for pp = 1:length(pSD.wk.LFP1(1,:))
        loglog(PSD.auxF,pSD.wk.LFP1(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     loglog(PSD.auxF,pSD.wk.LFP1(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     loglog(PSD.auxF,pSD.wk.LFP1(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.wk.LFP1(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.wk.LFP1(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.wk.LFP1(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Power.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('mPFC - WK')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    ylim([0.000001 0.1])
    xlim([0.5,100])
    
    %% CA11
    % Phase NREM
    subplot(3,3,4)
    for pp = 1:length(pSD.nrem.LFP2(1,:))
        loglog(PSD.auxF,pSD.nrem.LFP2(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     loglog(PSD.auxF,pSD.nrem.LFP2(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     loglog(PSD.auxF,pSD.nrem.LFP2(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.nrem.LFP2(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.nrem.LFP2(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.nrem.LFP2(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Power.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('CA1-1 - NREM')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    ylim([0.000001 0.1])
    xlim([0.5,100])
    
    % pSD REM
    subplot(3,3,5)
    for pp = 1:length(pSD.rem.LFP2(1,:))
        loglog(PSD.auxF,pSD.rem.LFP2(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     loglog(PSD.auxF,pSD.rem.LFP2(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     loglog(PSD.auxF,pSD.rem.LFP2(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.rem.LFP2(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.rem.LFP2(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.rem.LFP2(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Power.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('CA1-1 - REM')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    ylim([0.000001 0.1])
    xlim([0.5,100])
    
    % Phase WK
    subplot(3,3,6)
    for pp = 1:length(pSD.wk.LFP2(1,:))
        loglog(PSD.auxF,pSD.wk.LFP2(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     loglog(PSD.auxF,pSD.wk.LFP2(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     loglog(PSD.auxF,pSD.wk.LFP2(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.wk.LFP2(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.wk.LFP2(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.wk.LFP2(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Power.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('CA1-1 - WK')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    ylim([0.000001 0.1])
    xlim([0.5,100])
    
    
    %% CA12
    % Phase NREM
    subplot(3,3,7)
    for pp = 1:length(pSD.nrem.LFP3(1,:))
        loglog(PSD.auxF,pSD.nrem.LFP3(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     loglog(PSD.auxF,pSD.nrem.LFP3(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     loglog(PSD.auxF,pSD.nrem.LFP3(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.nrem.LFP3(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.nrem.LFP3(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.nrem.LFP3(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Power.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('CA1-2 - NREM')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    ylim([0.000001 0.1])
    xlim([0.5,100])
    
    % pSD REM
    subplot(3,3,8)
    for pp = 1:length(pSD.rem.LFP3(1,:))
        loglog(PSD.auxF,pSD.rem.LFP3(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     loglog(PSD.auxF,pSD.rem.LFP3(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     loglog(PSD.auxF,pSD.rem.LFP3(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.rem.LFP3(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.rem.LFP3(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.rem.LFP3(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Power.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('CA1-2 - REM')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    ylim([0.000001 0.1])
    xlim([0.5,100])
    
    % Phase WK
    subplot(3,3,9)
    for pp = 1:length(pSD.wk.LFP3(1,:))
        loglog(PSD.auxF,pSD.wk.LFP3(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     loglog(PSD.auxF,pSD.wk.LFP3(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     loglog(PSD.auxF,pSD.wk.LFP3(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.wk.LFP3(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.wk.LFP3(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     loglog(PSD.auxF,pSD.wk.LFP3(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Power.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('CA1-2 - WK')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    ylim([0.000001 0.1])
    xlim([0.5,100])
    
    counter_loop = 0;
    
    % Save PDF
    sgtitle(sprintf('PSD for Days - %s',names{ii,jj}))
    print('-fillpage',fullfile(directory,'PSD',sprintf('PSD for Days - %s',names{ii,jj})),'-dpdf','-r0',ff)
    
    close(ff)
    
end
end