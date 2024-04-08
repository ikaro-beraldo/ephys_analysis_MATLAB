%% Plot PSD for each DAY
function plot_psds_compare_days_function(directory,names)

% Pre-allocate
pSD.nrem.LFP1 = nan(10001,size(names,1),size(names,2));
pSD.rem.LFP1 = nan(10001,size(names,1),size(names,2));
pSD.wk.LFP1 = nan(10001,size(names,1),size(names,2));

pSD.nrem.LFP2 = nan(10001,size(names,1),size(names,2));
pSD.rem.LFP2 = nan(10001,size(names,1),size(names,2));
pSD.wk.LFP2 = nan(10001,size(names,1),size(names,2));

pSD.nrem.LFP3 = nan(10001,size(names,1),size(names,2));
pSD.rem.LFP3 = nan(10001,size(names,1),size(names,2));
pSD.wk.LFP3 = nan(10001,size(names,1),size(names,2));

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
        pSD.nrem.LFP1(:,ii,jj) = nanmean(PSD.LFP1(:,GMM_NREM_All_Sort==1),2);
        pSD.rem.LFP1(:,ii,jj) = nanmean(PSD.LFP1(:,GMM_REM_All_Sort==1),2);
        pSD.wk.LFP1(:,ii,jj) = nanmean(PSD.LFP1(:,GMM_WK_All_Sort==1),2);

        pSD.nrem.LFP2(:,ii,jj) = nanmean(PSD.LFP2(:,GMM_NREM_All_Sort==1),2);
        pSD.rem.LFP2(:,ii,jj) = nanmean(PSD.LFP2(:,GMM_REM_All_Sort==1),2);
        pSD.wk.LFP2(:,ii,jj) = nanmean(PSD.LFP2(:,GMM_WK_All_Sort==1),2);

        pSD.nrem.LFP3(:,ii,jj) = nanmean(PSD.LFP3(:,GMM_NREM_All_Sort==1),2);
        pSD.rem.LFP3(:,ii,jj) = nanmean(PSD.LFP3(:,GMM_REM_All_Sort==1),2);
        pSD.wk.LFP3(:,ii,jj) = nanmean(PSD.LFP3(:,GMM_WK_All_Sort==1),2);

        clear Cxy_all Call
    end
end

%% Plot
color(1,:) = [214 214 214];
color(2,:) = [255 204 153];
color(3,:) = [255 128 0];
color(4,:) = [255 51 51];
color(5,:) = [153 0 0];
color = color./255;

% Figure options
options.handle     = figure('PaperSize', [29.7 21]);
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
%options.color_area = [243 169 114]./255;    % Orange theme
%options.color_line = [236 112  22]./255;
options.alpha      = 0.5;
options.line_width = 2;
options.error      = 'sem';
options.x_axis = PSD.auxF;

for jj = 1:size(names,2)
    % define the color plot
    options.color_area = color(jj,:);    % Blue theme
    options.color_line = color(jj,:);

    options.sp = subplot(3,3,1);
    hold on
    plot_areaerrorbar(pSD.nrem.LFP1(:,:,jj)',options)
    options.sp = subplot(3,3,2);
    hold on
    plot_areaerrorbar(pSD.rem.LFP1(:,:,jj)',options)
    options.sp = subplot(3,3,3);
    hold on
    plot_areaerrorbar(pSD.wk.LFP1(:,:,jj)',options)

    options.sp = subplot(3,3,4);
    hold on
    plot_areaerrorbar(pSD.nrem.LFP2(:,:,jj)',options)
    options.sp = subplot(3,3,5);
    hold on
    plot_areaerrorbar(pSD.rem.LFP2(:,:,jj)',options);
    options.sp = subplot(3,3,6);
    hold on
    plot_areaerrorbar(pSD.wk.LFP2(:,:,jj)',options)

    options.sp = subplot(3,3,7);
    hold on
    plot_areaerrorbar(pSD.nrem.LFP3(:,:,jj)',options)
    options.sp = subplot(3,3,8);
    hold on
    plot_areaerrorbar(pSD.rem.LFP3(:,:,jj)',options)
    options.sp = subplot(3,3,9);
    hold on
    plot_areaerrorbar(pSD.wk.LFP3(:,:,jj)',options)
end

sgtitle(sprintf('PSD for Days - %s',names{ii,jj}))
print('-fillpage',fullfile(directory,'PSD','PSD Days comparison'),'-dpdf','-vector','-r0',options.handle)

end