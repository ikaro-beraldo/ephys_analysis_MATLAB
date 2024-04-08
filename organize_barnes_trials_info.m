ID = str2num(task_parameters.ID);
day = str2num(task_parameters.Day);
trial = str2num(task_parameters.Trial);

% Get the values
latency = task_parameters.Latency;
p_err = task_parameters.P_error;
s_err = task_parameters.S_error;
strategy = task_parameters.Strategy;
distance = task_parameters.Distance;
speed = task_parameters.Av_speed;
time_on_target = task_parameters.Time_on_target;
time_on_each_quadrant = task_parameters.Time_each_quadrant;

% Pre-allocate
trial_latency = nan(5,16);
trial_p_err = nan(5,16);
trial_s_err = nan(5,16);
trial_strategy = nan(5,16);
trial_distance = nan(5,16);
trial_speed = nan(5,16);
trial_time_on_target = nan(5,16);
trial_time_on_each_quadrant = nan(5,16);

u_ID = unique(ID);
u_day = unique(day);
u_trial = unique(trial);

% Organize the results by trial
for i = 1:length(u_ID)
    for j = 1:length(u_day)
        a = (j-1)*4 + 1;
        for z = 1:length(u_trial)
            b = z-1;
            idx = a+b;

            idx_v = ID == u_ID(i) & day == u_day(j) & trial == u_trial(z);

            if sum(idx_v) == 0
                continue
            end

            trial_latency(i,idx) = latency(idx_v);
            trial_p_err(i,idx) = p_err(idx_v);
            trial_s_err(i,idx) = s_err(idx_v);
            trial_strategy(i,idx) = strategy(idx_v,2);
            trial_distance(i,idx) = distance(idx_v);
            trial_speed(i,idx) = speed(idx_v);
            trial_time_on_target(i,idx) = time_on_target(idx_v);
            trial_time_on_each_quadrant(i,idx) = time_on_each_quadrant(idx_v,1);

        end
    end
end

% Process strategy matrix
trial_strategy(trial_strategy == 97) = 3;
trial_strategy(trial_strategy == 101) = 2;
trial_strategy(trial_strategy == 112) = 1;

trial_strategy_final = nan(3,16);


for s = 1:16
    trial_strategy_final(1,s) = sum(trial_strategy(:,s) == 1);  % Spacial
    trial_strategy_final(2,s) = sum(trial_strategy(:,s) == 2);  % Serial
    trial_strategy_final(3,s) = sum(trial_strategy(:,s) == 3);  % Random


end

% Permute the values randomly
for i = 1:5
    perm_v = randperm(size(trial_latency,2));
    perm_trial_latency(i,:) = trial_latency(i,perm_v);
end

%% Organize the results by day

% Pre-allocate
day_latency = nan(5,4);
day_p_err = nan(5,4);
day_s_err = nan(5,4);
day_strategy = nan(3,4);
day_distance = nan(5,4);
day_speed = nan(5,4);
day_time_on_target = nan(5,4);
day_time_on_each_quadrant = nan(5,4);

u_ID = unique(ID);
u_day = unique(day);
u_trial = unique(trial);

% Organize the results by day
for i = 1:length(u_ID)
    for j = 1:length(u_day)

        idx_v = ID == u_ID(i) & day == u_day(j);

        if sum(idx_v) == 0
            continue
        end

        day_latency(i,j) = nanmean(latency(idx_v));
        day_p_err(i,j) = nansum(p_err(idx_v));
        day_s_err(i,j) = nansum(s_err(idx_v));
        day_distance(i,j) = nanmean(distance(idx_v));
        day_speed(i,j) = nanmean(speed(idx_v));
        day_time_on_target(i,j) = nanmean(time_on_target(idx_v));
        day_time_on_each_quadrant(i,j) = nanmean(time_on_each_quadrant(idx_v,1));
        day_strategy_direct(i,j) = sum(strategy(idx_v,3) == 'a');
        day_strategy_serial(i,j) = sum(strategy(idx_v,3) == 'r');
        day_strategy_random(i,j) = sum(strategy(idx_v,3) == 'n');


        % Only for strategy
        idx = day == u_day(j);
        day_strategy(1,j) = sum(strategy(idx,3) == 'a');  % Third letter of 'spatial'
        day_strategy(2,j) = sum(strategy(idx,3) == 'r');  % Third letter of 'serial'
        day_strategy(3,j)= sum(strategy(idx,3) == 'n');  % Third letter of 'random'

    end

end

