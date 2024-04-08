%% Plot PSD for each DAY
function plot_PSD_compare_states_first_day_function(directory,names,channel_selection)

% Pre-allocate
mpfc.nrem = nan(10001,size(names,1),size(names,2));
mpfc.rem = nan(10001,size(names,1),size(names,2));
mpfc.wk = nan(10001,size(names,1),size(names,2));

ca1.nrem = nan(10001,size(names,1),size(names,2));
ca1.rem = nan(10001,size(names,1),size(names,2));
ca1.wk = nan(10001,size(names,1),size(names,2));

counter_loop = 0;

for ii = 1:size(names,1)
    for jj = 1:size(names,2)

        channel = channel_selection(ii);

        % Load file
        open_psd = fullfile(directory,'PSD',char(names{ii,jj}));
        load(open_psd,'PSD')

        if channel == 2
            CA1_PSD = PSD.LFP2;
        else
            CA1_PSD = PSD.LFP3;
        end
        
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
        mpfc.nrem(:,ii,jj) = nanmean(PSD.LFP1(:,GMM_NREM_All_Sort==1),2);
        mpfc.rem(:,ii,jj) = nanmean(PSD.LFP1(:,GMM_REM_All_Sort==1),2);
        mpfc.wk(:,ii,jj) = nanmean(PSD.LFP1(:,GMM_WK_All_Sort==1),2);

        ca1.nrem(:,ii,jj) = nanmean(CA1_PSD(:,GMM_NREM_All_Sort==1),2);
        ca1.rem(:,ii,jj) = nanmean(CA1_PSD(:,GMM_REM_All_Sort==1),2);
        ca1.wk(:,ii,jj) = nanmean(CA1_PSD(:,GMM_WK_All_Sort==1),2);

        clear CA1_PSD PSD.LFP1
    end
end

%% Plot
color(3,:)=[0.9290, 0.6940, 0.1250];    % WK   
color(1,:) =[0 0.4470 0.7410];          % NREM
color(2,:) =[0.3 0.3 0.3];              % REM

% Figure options
options.handle     = figure;
set(options.handle, 'Renderer', 'painters')
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
%options.color_area = [243 169 114]./255;    % Orange theme
%options.color_line = [236 112  22]./255;
options.alpha      = 0.5;
options.line_width = 2;
options.error      = 'sem';
options.x_axis = PSD.auxF(PSD.auxF<=50);
idx = find(PSD.auxF<=50);


for jj = 1:size(names,2)
    % define the color plot
    options.color_area = color(jj,:);    % Blue theme
    options.color_line = color(jj,:);

    % Plot phase
    
    options.sp = subplot(1,2,1);
    hold on
    options.color_area = color(1,:);  
    options.color_line = color(1,:);
    plot_areaerrorbar(mpfc.nrem(idx,:,jj)',options)

    options.color_area = color(2,:);
    options.color_line = color(2,:);
    plot_areaerrorbar(mpfc.rem(idx,:,jj)',options)

    options.color_area = color(3,:);
    options.color_line = color(3,:);
    plot_areaerrorbar(mpfc.wk(idx,:,jj)',options)

%     legend('NREM','REM','WK')
    title('PSD mPFC')
    set(gca, 'YScale', 'log')
    xlim([0 50])
%     set(gca, 'XScale', 'log')
%     legend boxoff
%     box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);

    % Plot spectral

    options.sp = subplot(1,2,2);
    set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
    hold on
    options.color_area = color(1,:);  
    options.color_line = color(1,:);
    plot_areaerrorbar(ca1.nrem(idx,:,jj)',options)

    options.color_area = color(2,:);
    options.color_line = color(2,:);
    plot_areaerrorbar(ca1.rem(idx,:,jj)',options)

    options.color_area = color(3,:);
    options.color_line = color(3,:);
    plot_areaerrorbar(ca1.wk(idx,:,jj)',options)

%     legend('NREM','REM','WK')
    title('PSD CA1')
    xlim([0 50])
%     legend boxoff
%     box off
    set(gca,'fontsize',18)
    set(gca,'Tickdir','out')
    set(gca,'Linewidth',1.5)
    set(gcf,'color',[1 1 1]);

end

sgtitle(sprintf('PSD for Days - %s',names{ii,jj}))
print('-bestfit',fullfile(directory,'PSD Day 1'),'-dpdf','-vector','-r0',options.handle)

end