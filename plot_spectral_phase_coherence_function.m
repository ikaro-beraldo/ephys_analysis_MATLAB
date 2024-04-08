%% Plot spectral and phase coherence throughout days

%% Plotting Spectral Coherence and Phase Coherence

function plot_spectral_phase_coherence_function(directory,names)

counter_loop = 0;
aux_Call = zeros(103,2201);

for ii = 1:size(names,1)
    for jj = 1:size(names,2)
        
        % Load files
        open_phase = fullfile(directory,'phase_coherence',char(names{ii,jj}));
        load(open_phase)
        open_spect = fullfile(directory,'spectral_coherence',char(names{ii,jj}));
        load(open_spect)
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
        phase.nrem(:,jj) = nanmean(Call(:,GMM_NREM_All_Sort==1),2);
        phase.rem(:,jj) = nanmean(Call(:,GMM_REM_All_Sort==1),2);
        phase.wk(:,jj) = nanmean(Call(:,GMM_WK_All_Sort==1),2);
        
        spectral.nrem(:,jj) = nanmean(Cxy_all(:,GMM_NREM_All_Sort==1),2);
        spectral.rem(:,jj) = nanmean(Cxy_all(:,GMM_REM_All_Sort==1),2);
        spectral.wk(:,jj) = nanmean(Cxy_all(:,GMM_WK_All_Sort==1),2);
        
        isequal(Call,aux_Call)
        aux_Call = Call;
        
        
        clear Cxy_all Call
    end
    
    counter_loop = counter_loop + 1;
    if counter_loop == 1
        colors = get(0, 'defaultAxesColorOrder');
        %         f=figure('PaperSize', [60 42],'visible','on');
        ff=figure('PaperSize', [29.7 21]);
        counter = 0;
    end
    
    % Spectral NREM
    subplot(2,3,1)
    for pp = 1:length(spectral.nrem(1,:))
        plot(F,spectral.nrem(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     plot(F,spectral.nrem(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     plot(F,spectral.nrem(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(F,spectral.nrem(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(F,spectral.nrem(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(F,spectral.nrem(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Spec. Coh.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('Spectral Coherence - NREM')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    xlim([0.5 30])
    ylim([0.1 1])
    
    % Spectral REM
    subplot(2,3,2)
    for pp = 1:length(spectral.rem(1,:))
        plot(F,spectral.rem(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     plot(F,spectral.rem(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     plot(F,spectral.rem(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(F,spectral.rem(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(F,spectral.rem(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(F,spectral.rem(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Spec. Coh.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('Spectral Coherence - REM')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    xlim([0.5 30])
    ylim([0.1 1])
    
    % Spectral WK
    subplot(2,3,3)
    for pp = 1:length(spectral.wk(1,:))
        plot(F,spectral.wk(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     plot(F,spectral.wk(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     plot(F,spectral.wk(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(F,spectral.wk(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(F,spectral.wk(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(F,spectral.wk(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Spec. Coh.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('Spectral Coherence - WK')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    xlim([0.5 30])
    ylim([0.1 1])
    
    % Phase NREM
    subplot(2,3,4)
    for pp = 1:length(phase.nrem(1,:))
        plot(F,phase.nrem(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     plot(f,phase.nrem(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     plot(f,phase.nrem(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(f,phase.nrem(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(f,phase.nrem(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(f,phase.nrem(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Spec. Coh.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('Phase Coherence - NREM')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    xlim([0.5 30])
    ylim([0.1 1])
    
    % phase REM
    subplot(2,3,5)
    for pp = 1:length(phase.rem(1,:))
        plot(F,phase.rem(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     plot(f,phase.rem(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     plot(f,phase.rem(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(f,phase.rem(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(f,phase.rem(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(f,phase.rem(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Spec. Coh.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('Phase Coherence - REM')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    xlim([0.5 30])
    ylim([0.1 1])
    
    % Phase WK
    subplot(2,3,6)
    for pp = 1:length(phase.wk(1,:))
        plot(F,phase.wk(:,pp),'Color',colors(pp,:)*0.8,'linewidth',1); xlim([0.5,100]);
        hold on
    end
%     plot(f,phase.wk(:,1),'Color',colors(1,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     hold on
%     plot(f,phase.wk(:,2),'Color',colors(2,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(f,phase.wk(:,3),'Color',colors(3,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(f,phase.wk(:,4),'Color',colors(4,:)*0.8,'linewidth',1); xlim([0.5,100]);
%     plot(f,phase.wk(:,5),'Color',colors(5,:)*0.8,'linewidth',1); xlim([0.5,100]);
    ylabel('Spec. Coh.')
    xlabel('Frequency (Hz)')
    legend('D1','D2','D3','D4','D5')
    title('Phase Coherence - WK')
    legend boxoff
    box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);
    xlim([0.5 30])
    ylim([0.1 1])
    
    counter_loop = 0;
    
    % Save PDF
    sgtitle(sprintf('Spectral and Phase Coherence - %s',names{ii,jj}))
    print('-fillpage',fullfile(directory,'spectral_coherence',sprintf('Spectral and Phase Coherence - %s',names{ii,jj})),'-dpdf','-r0',ff)
    
end
end