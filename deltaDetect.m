function delta_parameters = deltaDetect(data_windowing, fs, sd_threshold, plot_fig)
% This function uses wavelet decomposition to obtain delta frequency 0.4 - 4 Hz
% 
% input: Epoch data
%        fs sample frequency of data
% output : Calculated parameters of delta
% 
% Natalia Espinosa 25/04/2023
% Edited 21/07/2023
% ------------------------------------------------------------------------
count = 1;
count2 = 1;
count3 = 1;
step = 1/fs;

delta        = nan(size(data_windowing.data_w,1),size(data_windowing.data_w,2)); % Delta wave detected
delta_zscore = nan(size(data_windowing.data_w,1),size(data_windowing.data_w,2)); % z-score of delta
delta_diff   = nan(size(data_windowing.data_w,1),size(data_windowing.data_w,2)-1); % Derivate of delta

fig_index = 2;
subplot_index = 1;

for i = 1: size(data_windowing.data_w, 1) % i -> Number of trials  
    
    delta(i,:)        = wavelet(data_windowing.data_w(i,:), fs);
    delta_zscore(i,:) = zscore(delta(i,:));
    delta_diff(i,:)   = diff(delta_zscore(i,:))/step;

    indx_zerocross = find(delta_diff(i,2:end).* delta_diff(i,1:end-1)<0 )+1;

    %       For evoked delta
    %      if (type_detect == 'e')
    %         indx_zerocross = indx_zerocross((data_windowing.time_w(indx_zerocross)>=-0.2 & data_windowing.time_w(indx_zerocross)<4));
    %      end
               
    
    if (isempty(indx_zerocross)==0)
       % Find peak of delta wave
       [~,indx_peak]    = findpeaks(delta_zscore(i,:));
       [~,~,int_zcross] = intersect(indx_peak,indx_zerocross);
    
   
    % ===== Delta detection =====
    
    % Check delta events for each peak
        for jj = 1:length(int_zcross)
            flag = 0;
            if (indx_peak(jj)>1000 && indx_peak(jj)<6000)
                flag = 1;
                indx_onset_delta = 1;
                indx_peak_delta = indx_peak(1);
                indx_offset_delta = indx_zerocross(2);
            end
                   
            if ((int_zcross(jj) > 1) && (int_zcross(jj) < length(indx_zerocross)))
                flag =1;
                indx_peak_delta   = indx_zerocross(int_zcross(jj));     % Index vector of peak points
                indx_onset_delta  = indx_zerocross(int_zcross(jj) - 1); % Index vector of onset points
                indx_offset_delta = indx_zerocross(int_zcross(jj) + 1); % Index vector of offset points
            
            end
            if (flag==1)
                %plot(data_windowing.time_w(indx_onset_delta:indx_offset_delta), delta(i,indx_onset_delta:indx_offset_delta), 'LineWidth',1) 
                % ==== PARAMETERS =======
                %Amp_min     = 20; % mV
                Amp_min_sd = sd_threshold; % Standard Deviation
                Amplitude_onset   = abs(delta(i,indx_peak_delta) - delta(i,indx_onset_delta)); %[uV]
                Amplitude_offset   = abs(delta(i,indx_peak_delta) - delta(i,indx_offset_delta)); %[uV]
                Amplitude_sd_onset = abs(delta_zscore(i,indx_peak_delta) - delta_zscore(i,indx_onset_delta)); %[STD]
                Amplitude_sd_offset = abs(delta_zscore(i,indx_peak_delta) - delta_zscore(i,indx_offset_delta)); %[STD]
                time_peak   = data_windowing.time_w(indx_peak_delta);   % [seconds]
                time_onset  = data_windowing.time_w(indx_onset_delta);  % [seconds]
                time_offset = data_windowing.time_w(indx_offset_delta); % [seconds]
                Duration    = time_offset - time_onset;  % [seconds]
                latency = (time_peak - 0.08)*1000; % Latency of delta wave [ms] 

                % Slope onset
                Px = [time_onset; time_peak]; % Poinst of polygon in x
                Py = [delta(i,indx_onset_delta); delta(i,indx_peak_delta)]; % Poinst of polygon in y
                [coef_asscendet] = polyfit(Px,Py,1);
                Slope_onset      = coef_asscendet(1);
                
                % Slope offset
                Px_dec = [time_peak; time_offset]; % Poinst of polygon in x
                Py_dec = [delta(i,indx_peak_delta); delta(i,indx_offset_delta) ]; % Poinst of polygon in y
                [coef_desendent] = polyfit(Px_dec,Py_dec,1);
                Slope_offset     = coef_desendent(1);
                
                % Area under curve
                x_lim = data_windowing.time_w(indx_onset_delta:indx_offset_delta);
                Area = trapz(x_lim, delta(i,indx_onset_delta:indx_offset_delta)); % Area under curve [uV^2]
                Energy = trapz(x_lim, (delta(i,indx_onset_delta:indx_offset_delta)).^2); % Signal Energy
                
                
                % Well defined delta wave Amplitude minimum of 100 uV and time
                % onset minimum 50ms before event
                %                 if ((Amplitude > Amp_min) && (data_windowing.time_w(indx_onset_delta) > -0.05))
                if (Amplitude_sd_onset > Amp_min_sd) && Amplitude_sd_offset > Amp_min_sd
%                                      
                    % All delta detected
                    delta_parameters.delta_detected(count,:) = [i, time_onset, time_peak, time_offset, indx_onset_delta, indx_peak_delta,indx_offset_delta];
                    delta_parameters.parameters(count,:)   = [i, Amplitude_onset, Duration, Slope_onset, abs(Slope_offset), Area, Energy, latency];
                    count = count + 1;
                end
                                       
            end
            
            
        end
        
    end
   
% Plot Delta detected delta waves
% If it is necessary to plot
   if (plot_fig==1)
       % Generate figures with n subplots each 
       if (subplot_index > 12)
           subplot_index = 1;
           fig_index  = fig_index + 1;
           pause(0.5)
           close all
        end
       
        figure(fig_index)
        subplot(4,3,subplot_index)     
          
        %plot(data_windowing.time_w,data_windowing.data_w(i,:),'Color',[171 178 185]./255 )
        plot(data_windowing.time_w,data_windowing.data_w(i,:),'Color',[171 178 185]./255 )
        title(['trial: ' num2str(i)])
        xlim([0 data_windowing.time_w(end)])
        ylim([-250 250])
        hold on
        %Block delta idx
        delta_idx = find(delta_parameters.delta_detected(:,1) == i);
        for j = 1:length(delta_idx)
            onset  = find(data_windowing.time_w == delta_parameters.delta_detected(delta_idx(j),2));
            offset = find(data_windowing.time_w == delta_parameters.delta_detected(delta_idx(j),4));
            plot(data_windowing.time_w(onset:offset), data_windowing.data_w(i,onset:offset), 'LineWidth',1)        
        end  
        subplot_index = subplot_index + 1;
        
   end
  
delta_parameters.delta =  delta;
end
    prior_trial = 0;
    flag = 0;
        
    % First and second delta
    for m = 1:size(delta_parameters.delta_detected,1)
        trial = delta_parameters.delta_detected(m,1);
        if (trial ~= prior_trial && delta_parameters.delta_detected(m,6)<4000)
            flag = 1;
            delta_parameters.first_delta(count2,:) = delta_parameters.delta_detected(m,:);
            delta_parameters.parm_first_delta(count2,:) = delta_parameters.parameters(m,:);
            count2 = count2 + 1;
        end

        
        if (trial == prior_trial && flag == 1 && delta_parameters.delta_detected(m,6)>4000 && delta_parameters.delta_detected(m,6)<10000 )
            delta_parameters.second_delta(count3,:) = delta_parameters.delta_detected(m,:);
            delta_parameters.parm_second_delta(count3,:)  = delta_parameters.parameters(m,:);
            flag = 0;
            count3 = count3 + 1;
        end
        prior_trial = trial;
    end
end
 
