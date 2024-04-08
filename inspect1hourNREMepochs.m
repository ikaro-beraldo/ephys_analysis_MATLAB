% Visual inspection of the first hour
clear all

% directory = str(base directory name)
directory = 'E:\Barnes Maze - Mestrad\dados matlab\blocked_data';
% names2 = cell(matrix with each subdirectory separatelly
% names2 = {'B1_D1','B1_D2','B1_D3','B1_D4','B1_D5';'B2_D1','B2_D2','B2_D3','B2_D4','B2_D5';...
%     'B3_D1','B3_D2','B3_D3','B3_D4','B3_D5';...
%     'B4_D1','B4_D2','B4_D3','B4_D4','B4_D5';'B6_D1','B6_D2','B6_D3','B6_D4','B6_D5'};

names2 = {'B2_D1';...
    'B3_D1';...
    };

time_limit2 = (2 * 3600) / 10;
time_limit1 = (1 * 3600) / 10;

for ii = 1:size(names2,1)
    for jj = 1:size(names2,2)

        disp(names2{ii,jj})

        % Load GMM.All_Sort
        %% 0 - Get the classification results
        open_class = fullfile(directory,char(names2{ii,jj}),'GMM_Classification.mat');
        backup_file = fullfile(directory,char(names2{ii,jj}),'GMM_Classification_BKP.mat');
        load(open_class,'GMM','GMM_NREM_All_Sort','GMM_REM_All_Sort','GMM_WK_All_Sort')


        if ~isfile(backup_file)
            % Create a backup file for the classification if it doesn't exist
            % yet
            copyfile(open_class,backup_file)
        end

        % Fix the classification
        NREM = find(GMM.All_Sort == 2);
        REM = find(GMM.All_Sort == 1);
        WK = find(GMM.All_Sort == 3);

        % Filter by time
        NREM_1 = NREM(NREM <= time_limit1);
        REM_1 = REM(REM <= time_limit1);
        WK_1 = WK(WK <= time_limit1);

        % Filter by time
        NREM_2 = NREM(NREM <= time_limit2);
        REM_2 = REM(REM <= time_limit2);
        WK_2 = WK(WK <= time_limit2);
        %% Load LFP and Accel

%         open_lfp = fullfile(directory,char(names2{ii,jj}),'blocked_data.mat');
%         load(open_lfp,'Accel','LFP1','LFP3')
% 
%         for ee = 1:length(NREM)
%             idx = NREM(ee);
%             % Plot the data and Classification
%             subplot(3,1,1)
%             plot(LFP1(idx,:))
%             xlim([0 10000])
%             ylim([-1000 1000])
%             subplot(3,1,2)
%             plot(LFP3(idx,:))
%              xlim([0 10000])
%             ylim([-1000 1000])
%             subplot(3,1,3)
%             plot(Accel(idx,:))
%             xlim([0 10000])
%             ylim([0 0.05])
% 
%             % Get the answer to the classification (if it is right or
%             % wrong)
%             [~,~,button]=ginput(1);
%             switch button
%                 case 48 % Zero --> it is not correct
%                     GMM.All_Sort(NREM(ee)) = 3;  % Change it to awake
%                 case 51 % 3 --> AWAKE
%                     GMM.All_Sort(NREM(ee)) = 3;  % Change it to awake
%                 case 32 % Space --> Correct, go to the next one
%                 case 49 % 1 -- REM
%                     GMM.All_Sort(NREM(ee)) = 1;
%                 case 50 % 1 -- NREM
%                     GMM.All_Sort(NREM(ee)) = 2;
%             end
% 
%         end

        %% Load the BAND Info

        open_bands = fullfile(directory,'PSD',[char(names2{ii,jj}),'.mat']);
        load(open_bands,'BANDS')
        
        close all

        hold on
        title(char(names2{ii,jj}))
        plot(WK_2,BANDS.Gamma.LFP1(WK_2),'.','Color','yellow')
        plot(NREM_2,BANDS.Gamma.LFP1(NREM_2),'.','Color','blue')
        plot(REM_2,BANDS.Gamma.LFP1(REM_2),'.','Color',[.5 .5 .5])

        % Get the Y threshold
        [x,y]=ginput(1);

        % NREM above Y is awake and REM below Y is awake
        aux_WK = NREM_1(BANDS.Gamma.LFP1(NREM_1) > y);
        aux_WK = [aux_WK; REM_1(BANDS.Gamma.LFP1(REM_1) < y)];

        % Update GMM.All_Sort
        GMM.All_Sort(aux_WK) = 3;

        % Update the classification vectors
        GMM_REM_All_Sort = zeros(length(GMM_REM_All_Sort), 1);
        GMM_NREM_All_Sort = zeros(length(GMM_NREM_All_Sort), 1);
        GMM_WK_All_Sort = zeros(length(GMM_WK_All_Sort), 1);

        GMM_REM_All_Sort(GMM.All_Sort == 1) = 1;
        GMM_NREM_All_Sort(GMM.All_Sort == 2) = 1;
        GMM_WK_All_Sort(GMM.All_Sort == 3) = 1;

        % Save it
        save(open_class,"GMM","GMM_REM_All_Sort","GMM_NREM_All_Sort","GMM_WK_All_Sort",'-append')
    end
end