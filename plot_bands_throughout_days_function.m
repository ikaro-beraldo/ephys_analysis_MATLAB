%% Plot Bands throughout the days

%% MPFC
function plot_bands_throughout_days_function(directory,names,n_days)
% 3 figures for each animal (5 subplots: 1 for each band along the 5 days)

counter_loop = 0;

for ii = 1:length(names)
    
    counter_loop = counter_loop + 1;
    if counter_loop == 1
        colors = get(0, 'defaultAxesColorOrder');
        %         f=figure('PaperSize', [60 42],'visible','on');
        f=figure();
        counter = 0;
    end
    
    open_psd = fullfile(directory,'PSD',char(names{ii}));
    load(open_psd,'BANDS')
    open_class = fullfile(directory,char(names{ii}),'GMM_Classification.mat');
    load(open_class,'GMM','GMM_NREM_All_Sort','GMM_REM_All_Sort','GMM_WK_All_Sort')
    
    GMM_NREM_All_Sort(:) = false;
    GMM_NREM_All_Sort(GMM.All_Sort == 2) = true;
    GMM_REM_All_Sort(:) = false;
    GMM_REM_All_Sort(GMM.All_Sort == 1) = true;
    GMM_WK_All_Sort(:) = false;
    GMM_WK_All_Sort(GMM.All_Sort == 3) = true;
    
    e_start = 1 + counter;
    e_end = length(BANDS.Delta.LFP1) + counter;
    epochs_nrem = find(GMM_NREM_All_Sort==1)+counter;
    epochs_rem = find(GMM_REM_All_Sort==1)+counter;
    epochs_wk = find(GMM_WK_All_Sort==1)+counter;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       MPFC
    % Delta
    subplot(5,1,1)
    plot(epochs_nrem,BANDS.Delta.LFP1(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Delta.LFP1(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Delta.LFP1(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Delta')
    %ylim([0.000001 0.1])
    
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    % Theta
    subplot(5,1,2)
    plot(epochs_nrem,BANDS.Theta.LFP1(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Theta.LFP1(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Theta.LFP1(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Theta')
    ylabel('Power','FontSize',16)
    %ylim([0.000001 0.1])
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    % Theta 2
    subplot(5,1,3)
    plot(epochs_nrem,BANDS.Theta2.LFP1(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Theta2.LFP1(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Theta2.LFP1(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Theta 2')
    %ylim([0.000001 0.1])
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    
    % Theta 1
    subplot(5,1,4)
    plot(epochs_nrem,BANDS.Theta1.LFP1(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Theta1.LFP1(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Theta1.LFP1(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Theta 1')
    %ylim([0.000001 0.1])
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    
    % Gamma
    subplot(5,1,5)
    plot(epochs_nrem,BANDS.Gamma.LFP1(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Gamma.LFP1(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Gamma.LFP1(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Gamma')
    %ylim([0.000001 0.1])
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    counter = counter + length(BANDS.Delta.LFP1);
    if counter_loop == n_days
        counter_loop = 0;
        sgtitle(sprintf('Bands - mPFC %s',names{ii}))
        print('-fillpage',fullfile(directory,sprintf('Bands %s mPFC',names{ii})),'-dpdf','-r0',f)
    end
    
end

%% CA1 1

% 3 figures for each animal (5 subplots: 1 for each band along the 5 days)


counter_loop = 0;

for ii = 1:length(names)
    
    counter_loop = counter_loop + 1;
    if counter_loop == 1
        colors = get(0, 'defaultAxesColorOrder');
        %         f=figure('PaperSize', [60 42],'visible','on');
        f=figure();
        counter = 0;
    end
    
    open_psd = fullfile(directory,'PSD',char(names{ii}));
    load(open_psd,'BANDS')
    open_class = fullfile(directory,char(names{ii}),'GMM_Classification.mat');
    load(open_class,'GMM','GMM_NREM_All_Sort','GMM_REM_All_Sort','GMM_WK_All_Sort')
    
    GMM_NREM_All_Sort(:) = false;
    GMM_NREM_All_Sort(GMM.All_Sort == 2) = true;
    GMM_REM_All_Sort(:) = false;
    GMM_REM_All_Sort(GMM.All_Sort == 1) = true;
    GMM_WK_All_Sort(:) = false;
    GMM_WK_All_Sort(GMM.All_Sort == 3) = true;
    
    e_start = 1 + counter;
    e_end = length(BANDS.Delta.LFP2) + counter;
    epochs_nrem = find(GMM_NREM_All_Sort==1)+counter;
    epochs_rem = find(GMM_REM_All_Sort==1)+counter;
    epochs_wk = find(GMM_WK_All_Sort==1)+counter;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       MPFC
    % Delta
    subplot(5,1,1)
    plot(epochs_nrem,BANDS.Delta.LFP2(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Delta.LFP2(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Delta.LFP2(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Delta')
    %ylim([0.000001 0.1])
    
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    % Theta
    subplot(5,1,2)
    plot(epochs_nrem,BANDS.Theta.LFP2(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Theta.LFP2(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Theta.LFP2(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Theta')
    ylabel('Power','FontSize',16)
    %ylim([0.000001 0.1])
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    % Theta 2
    subplot(5,1,3)
    plot(epochs_nrem,BANDS.Theta2.LFP2(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Theta2.LFP2(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Theta2.LFP2(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Theta 2')
    %ylim([0.000001 0.1])
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    
    % Theta 1
    subplot(5,1,4)
    plot(epochs_nrem,BANDS.Theta1.LFP2(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Theta1.LFP2(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Theta1.LFP2(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Theta 1')
    %ylim([0.000001 0.1])
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    
    % Gamma
    subplot(5,1,5)
    plot(epochs_nrem,BANDS.Gamma.LFP2(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Gamma.LFP2(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Gamma.LFP2(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Gamma')
    %ylim([0.000001 0.1])
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    counter = counter + length(BANDS.Delta.LFP2);
    if counter_loop == n_days
        counter_loop = 0;
        sgtitle(sprintf('Bands - CA1-1 %s',names{ii}))
        print('-fillpage',fullfile(directory,sprintf('Bands %s CA11',names{ii})),'-dpdf','-r0',f)
    end
    
    
end

%% CA1 2

% 3 figures for each animal (5 subplots: 1 for each band along the 5 days)

counter_loop = 0;

for ii = 1:length(names)
    
    counter_loop = counter_loop + 1;
    if counter_loop == 1
        colors = get(0, 'defaultAxesColorOrder');
        %         f=figure('PaperSize', [60 42],'visible','on');
        f=figure();
        counter = 0;
    end
    
    open_psd = fullfile(directory,'PSD',char(names{ii}));
    load(open_psd,'BANDS')
    open_class = fullfile(directory,char(names{ii}),'GMM_Classification.mat');
    load(open_class,'GMM','GMM_NREM_All_Sort','GMM_REM_All_Sort','GMM_WK_All_Sort')
    
    GMM_NREM_All_Sort(:) = false;
    GMM_NREM_All_Sort(GMM.All_Sort == 2) = true;
    GMM_REM_All_Sort(:) = false;
    GMM_REM_All_Sort(GMM.All_Sort == 1) = true;
    GMM_WK_All_Sort(:) = false;
    GMM_WK_All_Sort(GMM.All_Sort == 3) = true;
    
    e_start = 1 + counter;
    e_end = length(BANDS.Delta.LFP3) + counter;
    epochs_nrem = find(GMM_NREM_All_Sort==1)+counter;
    epochs_rem = find(GMM_REM_All_Sort==1)+counter;
    epochs_wk = find(GMM_WK_All_Sort==1)+counter;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       MPFC
    % Delta
    subplot(5,1,1)
    plot(epochs_nrem,BANDS.Delta.LFP3(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Delta.LFP3(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Delta.LFP3(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Delta')
    %ylim([0.000001 0.1])
    
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    % Theta
    subplot(5,1,2)
    plot(epochs_nrem,BANDS.Theta.LFP3(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Theta.LFP3(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Theta.LFP3(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Theta')
    ylabel('Power','FontSize',16)
    %ylim([0.000001 0.1])
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    % Theta 2
    subplot(5,1,3)
    plot(epochs_nrem,BANDS.Theta2.LFP3(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Theta2.LFP3(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Theta2.LFP3(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Theta 2')
    %ylim([0.000001 0.1])
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    
    % Theta 1
    subplot(5,1,4)
    plot(epochs_nrem,BANDS.Theta1.LFP3(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Theta1.LFP3(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Theta1.LFP3(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Theta 1')
    %ylim([0.000001 0.1])
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    
    % Gamma
    subplot(5,1,5)
    plot(epochs_nrem,BANDS.Gamma.LFP3(GMM_NREM_All_Sort==1),'.','MarkerSize',3,'Color',colors(1,:)); xlim([1,e_end]);
    hold on
    plot(epochs_rem,BANDS.Gamma.LFP3(GMM_REM_All_Sort==1),'.','MarkerSize',3,'Color',[.3 .3 .3]); xlim([1,e_end]);
    plot(epochs_wk,BANDS.Gamma.LFP3(GMM_WK_All_Sort==1),'.','MarkerSize',3,'Color',colors(3,:)); xlim([1,e_end]);
    xline(e_start)
    
    %     legend('NREM','REM','WK','FontSize',14)
    %     legend boxoff
    title('Gamma')
    %ylim([0.000001 0.1])
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('Power','FontSize',16)
    
    set(gcf,'color',[1 1 1]);
    box off
    set(gca,'fontsize',14)
    set(gca,'Tickdir','out')
    set(gca,'fontname','helvetica')
    set(gca,'Linewidth',1.5)
    
    counter = counter + length(BANDS.Delta.LFP3);
    if counter_loop == n_days
        counter_loop = 0;
        sgtitle(sprintf('Bands - CA1-2 %s',names{ii}))
        print('-fillpage',fullfile(directory,sprintf('Bands %s CA12',names{ii})),'-dpdf','-r0',f)
    end
    
    
end
close all
end