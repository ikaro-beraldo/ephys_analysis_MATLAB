clear all
for animal = [1 2 3 4 6]
    for day = [1 2 3 4 5]
        name = sprintf('B%d_D%d',animal,day);

        mi_delta = fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\MI_delta',name);
        classification = fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data',name,'GMM_Classification');
        data = fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data',name,'blocked_data');
        ripple = fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\ripple',name);

        load(mi_delta,'MI_all')
        load(classification,'GMM')
        if animal == 4
            load(data,'LFP2')
            LFP3 = LFP2;
            clear LFP2
        else
            load(data,'LFP3')
        end
        load(ripple,'ripple_accepted')

        %% Compare MI in epochs with Ripple and without
        %         NREM = find(GMM.All_Sort==2);
        %
        %         ep_ripple = unique(ripple_accepted(:,1));
        %         ep_non_ripple = NREM(find(~ismember(NREM,ep_ripple)));
        %
        %         mi_90_140 = sum(MI_all.LFP3(4:6,:));
        %
        %         histogram(mi_90_140(ep_non_ripple),'Normalization','probability')
        %         hold on
        %         histogram(mi_90_140(ep_ripple),'Normalization','probability')
        %
        %         mean(mi_90_140(ep_non_ripple))
        %         mean(mi_90_140(ep_ripple))
        %
        %         [~,~,ix] = unique(ripple_accepted(:,1));
        %         epoch_counts = accumarray(ix,1);
        %
        %         plot(ep_ripple,epoch_counts)
        %
        %         [R,p] = corr(epoch_counts,mi_90_140(ep_ripple)');
        %         figure
        %         scatter(epoch_counts,mi_90_140(ep_ripple))

        %
        % [~,idx_max] = max(MI_all.LFP3(4,NREM));
        % [~,idx_min] = min(MI_all.LFP3(4,NREM));
        %
        % min_lfp_filt = eegfilt2(LFP3(NREM(idx_min),:),1000,150,170);
        % max_lfp_filt = eegfilt2(LFP3(NREM(idx_max),:),1000,150,170);
        %
        % min_lfp = eegfilt2(LFP3(NREM(idx_min),:),1000,1,4);
        % max_lfp = eegfilt2(LFP3(NREM(idx_max),:),1000,1,4);
        %
        % min_lfp_n = LFP3(NREM(idx_min),:);
        % max_lfp_n = LFP3(NREM(idx_max),:);
        %
        % subplot(4,1,1)
        % plot(min_lfp_n)
        % subplot(4,1,2)
        % plot(min_lfp)
        % hold on
        % plot(min_lfp_filt)
        % subplot(4,1,3)
        % plot(max_lfp_n)
        % subplot(4,1,4)
        % plot(max_lfp)
        % hold on
        % plot(max_lfp_filt)
        %
        % % plot possible ripples
        % rps = find(ripple_accepted(:,1) == NREM(idx_max));
        % subplot(4,1,4)
        % for i = 1:length(rps)
        %    idx = NREM(idx_max);
        %    a = ripple_accepted(rps(i),2);
        %    b = ripple_accepted(rps(i),4);
        %
        %    plot(a:b,max_lfp_filt(a:b))
        % end
        %
        % figure
        % plot(angle(hilbert(max_lfp)))
        % figure
        % plot(abs(hilbert(max_lfp)))


        %%
        idx_aux = 0;
        rp_delta_phase = nan(size(ripple_accepted,1),1);
        total = size(ripple_accepted,1);
        fprintf('Total: %d\n',total)
        for i = 1:size(ripple_accepted,1)
            fprintf('%d.',i)

            idx = ripple_accepted(i,1);
            a = ripple_accepted(i,2);
            b = ripple_accepted(i,4);
            peak = ripple_accepted(i,3);

            if idx ~= idx_aux
                %rp_filt = eegfilt2(LFP3(idx,:),1000,100,220);
                del_filt = angle(hilbert(eegfilt2(LFP3(idx,:),1000,1,4)));
            end

            rp_delta_phase(i) = del_filt(peak);
            idx_aux = idx;
        end


        % Transform the radians to degree and compute entropy
        rp_delta_phase = rad2deg(rp_delta_phase) +180;
        entropy = -sum(rp_delta_phase.*log2(rp_delta_phase));

        f = figure;
        nbin = 18;
        histogram(rp_delta_phase(:,1),nbin,'Normalization','probability')
        % MI=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin)        

        signal = sin(-pi:0.01:pi);
        t = linspace(0,360,length(signal));
        hold on
        plot(t,signal.*0.05 + 0.05)

        title(sprintf('Name: %s; Entropy: %d',name,entropy))

        sav = fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\ripple_MI',name);
        save(sav,'rp_delta_phase')

    end
end

%%
clear all
for animal = [1 2 3 4 6]
    for day = [1 2 3 4 5]
        name = sprintf('B%d_D%d',animal,day);
        load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\ripple_MI',name),'rp_delta_phase');
        load(fullfile('E:\Barnes Maze - Mestrad\dados matlab\blocked_data\ripple',name),'ripple_accepted');

        % Found the ripples regarding the time limit
        limit = ripple_accepted(:,1) <= 180 * 60 /10;
        n_bins = 18;

        [N,edges,bin]=histcounts(rp_delta_phase(limit),n_bins);
%         rp_delta_phase_binned = circ_ang2rad((bin-1).*20);
        rad_edges = circ_ang2rad(edges);

        % Normalize so the sum is 1
        N_for_entropy = N./sum(N);

        entropy(animal,day) = -sum(N_for_entropy.*log2(N_for_entropy));

        % Convert the degree to radian
        rad_rp_delta_phase = circ_ang2rad(rp_delta_phase);
       
        if day == 1
            figure
        end

        h = subplot(1,5,day);
        % Polar Plot
        polarhistogram(rad_rp_delta_phase,rad_edges,'Normalization','probability','FaceAlpha',0)
        rlim([0 0.2])

        % Define the mean
        mean_angle = circ_mean(rad_rp_delta_phase);
        av_angle(animal,day) = mean_angle;
    
        hold on
        polarplot([mean_angle mean_angle], [0 max(N)], 'r', 'LineWidth', 2);
        hold off;

        %         f = figure;
        %         nbin = 18;
        %         histogram(rp_delta_phase(:,1),nbin,'Normalization','probability')
        %         % MI=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin)
        %
        %         signal = sin(-pi:0.01:pi);
        %         t = linspace(0,360,length(signal));
        %         hold on
        %         plot(t,signal.*0.05 + 0.05)
        %
        %         title(sprintf('Name: %s; Entropy: %d',name,entropy))
    end
end


% alpha = randn(60,1)*.4+pi/2;
%     figure
%     subplot(2,2,1)
%     circ_plot(alpha,'pretty','ro',true,'linewidth',2,'color','r'),
%     title('pretty plot style')
%     subplot(2,2,2)
%     circ_plot(alpha,'hist',[],20,true,true,'linewidth',2,'color','r')
%     title('hist plot style')
%     subplot(2,2,3)
%     circ_plot(alpha,[],'s')
%     title('non-fancy plot style')