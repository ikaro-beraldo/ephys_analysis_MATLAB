entropy.delta.nrem = [];
entropy.delta.rem = [];
entropy.delta.wk = [];

entropy.theta.nrem = [];
entropy.theta.rem = [];
entropy.theta.wk = [];

n_blocks = 3;
time_limit = 1 : 360 * 60/10

for ii = [1 2 3 4 6]
    for jj = 1:5
        name = sprintf('B%d_D%d',ii,jj);
        figure
        load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\comodulation\CA1',name))
        load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data',name,'GMM_Classification.mat'),'GMM')

        phase_limit_delta = find(parameters.PhaseFreqVector >= 0.5 & parameters.PhaseFreqVector <=4);
        phase_limit_theta = find(parameters.PhaseFreqVector > 4 & parameters.PhaseFreqVector <=10);

        Comodulogram_theta = permute(sum(Comodulogram_all(phase_limit_theta,:,:)),[2 3 1])';
        Comodulogram_delta = permute(sum(Comodulogram_all(phase_limit_delta,:,:)),[2 3 1])';

        spec_Entropy_delta = -sum(Comodulogram_delta.*log2(Comodulogram_delta),2);
        spec_Entropy_theta = -sum(Comodulogram_theta.*log2(Comodulogram_theta),2);

        REM = find(GMM.All_Sort==1);
        NREM = find(GMM.All_Sort==2);
        WK = find(GMM.All_Sort==3);

%         % Intersect the states inside the time limit
%         REM = intersect(REM, time_limit);
%         NREM = intersect(NREM, time_limit);
%         WK = intersect(WK, time_limit);
       
         % Define exactly the blocks of stable states
        [REM_idx, NREM_idx, WK_idx] = define_beg_end_for_n_blocks(GMM.All_Sort,n_blocks);

        % Intersect the states inside the time limit
        REM = intersect(REM_idx, time_limit);
        NREM = intersect(NREM_idx, time_limit);
        WK = intersect(WK_idx, time_limit);

        color(3,:)=[0.9290, 0.6940, 0.1250];    % WK
        color(1,:) =[0 0.4470 0.7410];          % NREM
        color(2,:) =[0.3 0.3 0.3];              % REM

%         Phase Delta
        subplot(2,1,1)
        hold on
        plot(NREM,spec_Entropy_delta(NREM),'.','Color',color(1,:))
        plot(REM,spec_Entropy_delta(REM),'.','Color',color(2,:))
        plot(WK,spec_Entropy_delta(WK),'.','Color',color(3,:))
        ylim([0 10])
        xlim([0 max([NREM; REM; WK])])
        title('Phase Coherence - Delta')
        ylabel('Coherence')



        set(gca,'Tickdir','out')
        set(gca,'Linewidth',1.5)
        set(gcf,'color',[1 1 1]);

        % Phase Theta
        subplot(2,1,2)
        hold on
        plot(NREM,spec_Entropy_theta(NREM),'.','Color',color(1,:))
        plot(REM,spec_Entropy_theta(REM),'.','Color',color(2,:))
        plot(WK,spec_Entropy_theta(WK),'.','Color',color(3,:))
        ylim([0 10])
        xlim([0 max([NREM; REM; WK])])
        title('Phase Coherence - Theta')
        ylabel('Coherence')



        set(gca,'Tickdir','out')
        set(gca,'Linewidth',1.5)
        set(gcf,'color',[1 1 1]);

        % Change the ii to the proper value
        if ii == 6
            ii = 5;
        end

        % Save values
        entropy.delta.nrem(ii,jj) = nanmean(spec_Entropy_delta(NREM));
        entropy.delta.rem(ii,jj) = nanmean(spec_Entropy_delta(REM));
        entropy.delta.wk(ii,jj) = nanmean(spec_Entropy_delta(WK));

        entropy.theta.nrem(ii,jj) = nanmean(spec_Entropy_theta(NREM));
        entropy.theta.rem(ii,jj) = nanmean(spec_Entropy_theta(REM));
        entropy.theta.wk(ii,jj) = nanmean(spec_Entropy_theta(WK));
        
        if ii == 5
            ii = 6;
        end
    end
end

function [REM_final, NREM_final, WK_final] = define_beg_end_for_n_blocks(All_Sort,n_blocks)

    % Get REM blocks
    REM = zeros(1,length(All_Sort));
    REM(All_Sort == 1) = 1;
    REM_idx = get_b_and_g(REM, n_blocks);

    % Get NREM blocks
    NREM = zeros(1,length(All_Sort));
    NREM(All_Sort == 2) = 1;
    NREM_idx = get_b_and_g(NREM, n_blocks);

    % Get WK blocks
    WK = zeros(1,length(All_Sort));
    WK(All_Sort == 3) = 1;
    WK_idx = get_b_and_g(WK, n_blocks);

    % Create the final NREM vector
    NREM_final = nan(sum((NREM_idx(:,2) - NREM_idx(:,1))+1),1);

    init = 1;
    for i = 1:size(NREM_idx,1)        
        n = (NREM_idx(i,2) - NREM_idx(i,1));
        NREM_final(init:n+init) = NREM_idx(i,1) : NREM_idx(i,2);
        init = init + n + 1;
    end

    % Create the final REM vector
    REM_final = nan(sum((REM_idx(:,2) - REM_idx(:,1))+1),1);

    init = 1;
    for i = 1:size(REM_idx,1)        
        n = (REM_idx(i,2) - REM_idx(i,1));
        REM_final(init:n+init) = REM_idx(i,1) : REM_idx(i,2);
        init = init + n + 1;
    end


    % Create the final WK vector
    WK_final = nan(sum((WK_idx(:,2) - WK_idx(:,1))+1),1);

    init = 1;
    for i = 1:size(WK_idx,1)        
        n = (WK_idx(i,2) - WK_idx(i,1));
        WK_final(init:n+init) = WK_idx(i,1) : WK_idx(i,2);
        init = init + n + 1;
    end
    
end

% Define the beginning and end of each block
function state_idx = get_b_and_g(state, n_blocks)
primary_detection_b = find(diff(state)>0)+1;
primary_detection_e = find(diff(state)<0);

if primary_detection_b(1)>primary_detection_e(1)
    primary_detection_e = primary_detection_e(2:end);
end

if primary_detection_b(end)>primary_detection_e(end)
    primary_detection_b = primary_detection_b(1:end-1);
end

selected_idx = find(((primary_detection_e - primary_detection_b)+1) >= n_blocks);
state_idx = [primary_detection_b(selected_idx)' primary_detection_e(selected_idx)'];
end
