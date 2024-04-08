%% CORRELATION PROB TRIAL STATISTICS
% Load the Barnes Maze Prob Trial statistics
load("E:\Barnes Maze - Mestrad\Resultados DLC - Barnes\ProbTrial\Final_results.mat")
prob_time_on_each_quadrant = task_parameters.("Time_each_quadrant");
time_on_target_min1 = task_parameters.('Time_on_target_minute');

result_without_barnes = results;
result_without_barnes = rmfield(result_without_barnes,'barnes');
result_without_barnes = rmfield(result_without_barnes,'sleep_parameters');

% Get every single results field
results_names = fieldnamesr(result_without_barnes,'prefix');
correlation_results_prob = cell(length(results_names),4);
counter = 1;
for jj = 1:length(results_names)
    % Get the parameters
    %     a = prob_time_on_each_quadrant(:,1);
    a = time_on_target_min1(:,1);
    eval(sprintf('b = %s;',results_names{jj}));

    % Take the only the first hour of each day
        b = permute(b,[3 2 1]);
        b1 = b(1,:,1); % D1
        b2 = b(1,:,2) ;


        nan_values = isnan(b1) | isnan(b2);

        if(sum(nan_values) >= 1)
            b1(nan_values) = [];
            b2(nan_values) = [];
            a(nan_values) = [];
        end


    % Take the test parameter and average its days, except the baseline
%     b1 = nanmean(b(1,:,1:2),3);
%     b2 = nanmean(b(2,:,1:2),3);

    % Take only the first day of the task
    b = b2 ./ b1;
    %       b = b2;


    [correlation_results_prob{counter,1},correlation_results_prob{counter,2}] = corrcoef(a,b);
    % Fix
    if length(correlation_results_prob{counter,1}) == 1 % If the correlation results yields only 1 R and 1 p value
        correlation_results_prob{counter,1} = correlation_results_prob{counter,1}(1);
        correlation_results_prob{counter,2} = correlation_results_prob{counter,2}(1);
    else
        correlation_results_prob{counter,1} = correlation_results_prob{counter,1}(2);
        correlation_results_prob{counter,2} = correlation_results_prob{counter,2}(2);
    end
    correlation_results_prob{counter,3} = 'prob_time_each_quadrant';
    correlation_results_prob{counter,4} = results_names{jj};

    counter = counter + 1;
end


%% Bonferroni Correction
% n_tests = length([correlation_results_prob{:,2}]);
% for ii = 1:n_tests
%     correlation_results_prob{ii,2} = correlation_results_prob{ii,2}*n_tests;
% end

%% Filter correlation results

clear filtered_correlation
alfa = 0.05;
filtered = find([correlation_results_prob{:,2}] < alfa & [correlation_results_prob{:,2}] ~= 0 & [correlation_results_prob{:,1}] ~= 1);

for ii = 1:length(filtered)
    filtered_correlation{ii,1} = correlation_results_prob{filtered(ii),1};
    filtered_correlation{ii,2} = correlation_results_prob{filtered(ii),2};
    filtered_correlation{ii,3} = correlation_results_prob{filtered(ii),3};
    filtered_correlation{ii,4} = correlation_results_prob{filtered(ii),4};
end

%% Organize the results in a cell

final_correlation = cell(length(filtered_correlation),2);

for i = 1:length(correlation_results_prob)
    result_n = correlation_results_prob{i,4};
    % Prepare the data
    eval(sprintf('aux = %s;',result_n))
    aux = permute(aux,[3 2 1]);
    b1 = aux(1,:,1); % D1
    b2 = aux(1,:,2);

    % Take the test parameter and average its days, except the baseline
    %     b1 = nanmean(aux(1,:,1:2),3);
    %     b2 = nanmean(aux(2,:,1:2),3);

    % Take the test parameter and average its days, except the baseline
    %b = nanmean(b(2:end,:));
    % Take only the first day of the task
    aux = b2 ./ b1;
    
    final_correlation{i,1} = result_n;
    final_correlation{i,2} = [time_on_target_min1(:,1) aux'];

end