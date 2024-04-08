function get_MI_delta_degree_amplitude_epoch_is_N(directory,names,time_limit,b_length)

counter_loop = 0;
results = struct;

% Time limit to consider in minutes
time_limit = time_limit * 60 / b_length;


% for jj = 1:size(names,2)
for jj = [1 2 3 4 5]
    MI_amp.LFP2.f_30_55.NREM = [];
    MI_amp.LFP3.f_90_140.NREM = [];

    for ii = 1:size(names,1)

        %% 0 - Get the classification results
        open_class = fullfile(directory,char(names{ii,jj}),'GMM_Classification.mat');
        load(open_class,'GMM','GMM_NREM_All_Sort','GMM_REM_All_Sort','GMM_WK_All_Sort')

        % Fix the classification
        NREM = find(GMM.All_Sort == 2);
        REM = find(GMM.All_Sort == 1);
        WK = find(GMM.All_Sort == 3);

        % Filter by time
        NREM = NREM(NREM <= time_limit);
        REM = REM(REM <= time_limit);
        WK = WK(WK <= time_limit);


        % Load files
        filename = fullfile(directory,'MI_delta',char(names{ii,jj}));
        load(filename,'MeanAmp_all')

        MeanAmp_all.LFP2 = gather(MeanAmp_all.LFP2);
        MeanAmp_all.LFP3 = gather(MeanAmp_all.LFP3);

        Pf = [5 12];
        % Amplitude
        Af = [30 55; 65 90; 90 140; 150 170; 180 200; 200 250];

        % LFP2 - 30 - 55
        MI_amp.LFP2.f_30_55.NREM = cat(1,MI_amp.LFP2.f_30_55.NREM,reshape(MeanAmp_all.LFP2(1,NREM,:),[],18));
        MI_amp.LFP2.f_30_55.REM(jj,ii,:) = mean(MeanAmp_all.LFP2(1,REM,:),2,'omitnan');
        MI_amp.LFP2.f_30_55.WK(jj,ii,:) = mean(MeanAmp_all.LFP2(1,WK,:),2,'omitnan');

        % LFP2 - 65 - 90
        MI_amp.LFP2.f_65_90.NREM(jj,ii,:) = mean(MeanAmp_all.LFP2(2,NREM,:),2,'omitnan');
        MI_amp.LFP2.f_65_90.REM(jj,ii,:) = mean(MeanAmp_all.LFP2(2,REM,:),2,'omitnan');
        MI_amp.LFP2.f_65_90.WK(jj,ii,:) = mean(MeanAmp_all.LFP2(2,WK,:),2,'omitnan');

        % LFP2 - 90 - 140
        MI_amp.LFP2.f_90_140.NREM(jj,ii,:) = mean(MeanAmp_all.LFP2(3,NREM,:),2,'omitnan');
        MI_amp.LFP2.f_90_140.REM(jj,ii,:) = mean(MeanAmp_all.LFP2(3,REM,:),2,'omitnan');
        MI_amp.LFP2.f_90_140.WK(jj,ii,:) = mean(MeanAmp_all.LFP2(3,WK,:),2,'omitnan');

        % LFP2 - 150 - 170
        MI_amp.LFP2.f_150_170.NREM(jj,ii,:) = mean(MeanAmp_all.LFP2(4,NREM,:),2,'omitnan');
        MI_amp.LFP2.f_150_170.REM(jj,ii,:) = mean(MeanAmp_all.LFP2(4,REM,:),2,'omitnan');
        MI_amp.LFP2.f_150_170.WK(jj,ii,:) = mean(MeanAmp_all.LFP2(4,WK,:),2,'omitnan');

        % LFP2 - 180 - 200
        MI_amp.LFP2.f_180_200.NREM(jj,ii,:) = mean(MeanAmp_all.LFP2(5,NREM,:),2,'omitnan');
        MI_amp.LFP2.f_180_200.REM(jj,ii,:) = mean(MeanAmp_all.LFP2(5,REM,:),2,'omitnan');
        MI_amp.LFP2.f_180_200.WK(jj,ii,:) = mean(MeanAmp_all.LFP2(5,WK,:),2,'omitnan');

        % LFP2 - 200 - 250
        MI_amp.LFP2.f_200_250.NREM(jj,ii,:) = mean(MeanAmp_all.LFP2(6,NREM,:),2,'omitnan');
        MI_amp.LFP2.f_200_250.REM(jj,ii,:) = mean(MeanAmp_all.LFP2(6,REM,:),2,'omitnan');
        MI_amp.LFP2.f_200_250.WK(jj,ii,:) = mean(MeanAmp_all.LFP2(6,WK,:),2,'omitnan');

        %% LFP3

        % LFP2 - 30 - 55
        MI_amp.LFP3.f_30_55.NREM(jj,ii,:) = mean(MeanAmp_all.LFP3(1,NREM,:),2,'omitnan');
        MI_amp.LFP3.f_30_55.REM(jj,ii,:) = mean(MeanAmp_all.LFP3(1,REM,:),2,'omitnan');
        MI_amp.LFP3.f_30_55.WK(jj,ii,:) = mean(MeanAmp_all.LFP3(1,WK,:),2,'omitnan');

        % LFP3 - 65 - 90
        MI_amp.LFP3.f_65_90.NREM(jj,ii,:) = mean(MeanAmp_all.LFP3(2,NREM,:),2,'omitnan');
        MI_amp.LFP3.f_65_90.REM(jj,ii,:) = mean(MeanAmp_all.LFP3(2,REM,:),2,'omitnan');
        MI_amp.LFP3.f_65_90.WK(jj,ii,:) = mean(MeanAmp_all.LFP3(2,WK,:),2,'omitnan');

        % LFP3 - 90 - 140
        MI_amp.LFP3.f_90_140.NREM = cat(1,MI_amp.LFP3.f_90_140.NREM, reshape(MeanAmp_all.LFP3(3,NREM,:),[],18));

        MI_amp.LFP3.f_90_140.REM(jj,ii,:) = mean(MeanAmp_all.LFP3(3,REM,:),2,'omitnan');
        MI_amp.LFP3.f_90_140.WK(jj,ii,:) = mean(MeanAmp_all.LFP3(3,WK,:),2,'omitnan');

        % LFP3 - 150 - 170
        MI_amp.LFP3.f_150_170.NREM(jj,ii,:) = mean(MeanAmp_all.LFP3(4,NREM,:),2,'omitnan');
        MI_amp.LFP3.f_150_170.REM(jj,ii,:) = mean(MeanAmp_all.LFP3(4,REM,:),2,'omitnan');
        MI_amp.LFP3.f_150_170.WK(jj,ii,:) = mean(MeanAmp_all.LFP3(4,WK,:),2,'omitnan');

        % LFP3 - 180 - 200
        MI_amp.LFP3.f_180_200.NREM(jj,ii,:) = mean(MeanAmp_all.LFP3(5,NREM,:),2,'omitnan');
        MI_amp.LFP3.f_180_200.REM(jj,ii,:) = mean(MeanAmp_all.LFP3(5,REM,:),2,'omitnan');
        MI_amp.LFP3.f_180_200.WK(jj,ii,:) = mean(MeanAmp_all.LFP3(5,WK,:),2,'omitnan');

        % LFP3 - 200 - 250
        MI_amp.LFP3.f_200_250.NREM(jj,ii,:) = mean(MeanAmp_all.LFP3(6,NREM,:),2,'omitnan');
        MI_amp.LFP3.f_200_250.REM(jj,ii,:) = mean(MeanAmp_all.LFP3(6,REM,:),2,'omitnan');
        MI_amp.LFP3.f_200_250.WK(jj,ii,:) = mean(MeanAmp_all.LFP3(6,WK,:),2,'omitnan');


    end

%% Normalize the amplitude
% MI_amp.LFP3.f_90_140.NREM = mean(MI_amp.LFP3.f_90_140.NREM);
MI_amp.LFP3.f_90_140.NREM = MI_amp.LFP3.f_90_140.NREM./max(MI_amp.LFP3.f_90_140.NREM,[],2);

    %% Plot
color(1,:) = [214 214 214];
color(2,:) = [255 204 153];
color(3,:) = [255 128 0];
color(4,:) = [255 51 51];
color(5,:) = [153 0 0];
color = color./255;

% Figure options
options.handle     = figure(2);
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
%options.color_area = [243 169 114]./255;    % Orange theme
%options.color_line = [236 112  22]./255;
options.alpha      = 0.5;
options.line_width = 2;
options.error      = 'sem';
options.x_axis = 10:20:720;
options.sp = subplot(1,1,1);

% define the color plot
options.color_area = color(jj,:);    % Blue theme
options.color_line = color(jj,:);

hold on
% plot(10:20:720,[mean(MI_amp.LFP3.f_90_140.NREM) mean(MI_amp.LFP3.f_90_140.NREM)],"Color",color(jj,:))
plot_areaerrorbar([MI_amp.LFP3.f_90_140.NREM MI_amp.LFP3.f_90_140.NREM],options)
xlim([0 720])
nbin = 18;
MeanAmp = mean(MI_amp.LFP3.f_90_140.NREM);
MI(jj)=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);

end
MI
end