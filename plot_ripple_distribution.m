k = 1;
for ii = 1:5
    for jj = 1:6

        teste = (results_distribution_bin.RDS.duration.ripple{ii,jj});
        idx = kmeans(teste,1);

        %         figure(1)
        %         subplot(5,6,k)
        %         histogram(teste)
        % %         xlim([0 60])
        %

        %         figure(2)
        subplot(5,6,k)
        hold on
        histogram(teste(idx==1))
        histogram(teste(idx==2))
        histogram(teste(idx==3))
                 xlim([0 0.2])
        ylim([0 150])
        k = k + 1;
        av1 = round(median(teste(idx==1)),2);
        av2 = round(median(teste(idx==2)),2);
        text(av1,100,string(av1))
        text(av2,100,string(av2))

        %         load()

    end
end


%%

k = 1;
for ii = 1:5

    teste = zscore(results_distribution.RDS.RMS.ripple_2{ii,1});
    idx = kmeans(teste,2);

    %         figure(1)
    %         subplot(5,6,k)
    %         histogram(teste)
    % %         xlim([0 60])
    %


    %         figure(2)
    subplot(1,5,k)
    hold on
    histogram(teste(idx==1))
    histogram(teste(idx==2))
    histogram(teste(idx==3))
    %         xlim([0 60])
    ylim([0 700])
    k = k + 1;
    av1 = round(median(teste(idx==1)),2);
    av2 = round(median(teste(idx==2)),2);
    text(av1,600,string(av1))
    text(av2,600,string(av2))

    %         load()

end


%% 

names2 = {'B1_D1','B1_D2','B1_D3','B1_D4','B1_D5';'B2_D1','B2_D2','B2_D3','B2_D4','B2_D5';...
    'B3_D1','B3_D2','B3_D3','B3_D4','B3_D5';...
    'B4_D1','B4_D2','B4_D3','B4_D4','B4_D5';'B6_D1','B6_D2','B6_D3','B6_D4','B6_D5'};

k = 1;
av = [];
stds = [];
for ii = 1:5

    dist = [];

    for jj = 1:5

        file_name = fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\ripple',names2{ii,jj});
        load(file_name,'RMS')

        dist = [dist RMS.ripple_3];
    end

    av(ii) = mean(dist);
    stds(ii) = std(dist);

end

k = 1;
for ii = 1:5
    for jj = 1:5

        file_name = fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\ripple',names2{ii,jj});
        load(file_name,'RMS')

        teste = (RMS.ripple_3-av(ii))/stds(ii);
        idx = kmeans(teste',2);

        figure(1)
        subplot(5,5,k)
        histogram(teste)
        xlim([-4 4])
        ylim([0 800])
        av1 = round(median(teste),3);
        text(av1,700,string(av1))


        figure(2)
        subplot(5,5,k)
        hold on
        histogram(teste(idx==1))
        histogram(teste(idx==2))
        histogram(teste(idx==3))
        %         xlim([0 60])
        ylim([0 700])
        k = k + 1;
        av1 = round(median(teste(idx==1)),2);
        av2 = round(median(teste(idx==2)),2);
        text(av1,600,string(av1))
        text(av2,600,string(av2))

        %         load()

    end
end


