function [curve_fit_results,curve_fit_anova_results,curve_fit_anova_rejected, curve_fit_results_names, ttest_results, filtered_ttest] = run_curve_fit(results, plot_bool, curve_fit_mode)

% Get every single results field
results_names = fieldnamesr(results,'prefix');

for ii = 1:length(results_names)
    cm = sprintf('matrix = %s;',results_names{ii});
    % execute command
    eval(cm)

    % If the matrix is 3d create the curve fit
    if ndims(matrix) == 3
        % Permute the matrix so the days become hours and hours turn into days
        % (poetic)
        matrix = permute(matrix, [3 2 1]);
        % Get the curve fit
        result_curve_fit = curve_fit_ikaro(matrix, curve_fit_mode);

        if plot_bool
            % Plot the curve fit
            plot_curve_fit_multiple_days(matrix, results_names{ii}, curve_fit_mode)
        end

    elseif ndims(matrix) == 2   % if not, just ignore it
        continue
    end

    cm = sprintf('curve_fit_%s = result_curve_fit;', results_names{ii});
    eval(cm)
end

% Run anova to the curve_fit results
[curve_fit_anova_results,curve_fit_anova_rejected] = calculate_anova(curve_fit_results);

% Run a paired ttest
[ttest_results, filtered_ttest] = run_ttest(curve_fit_results);

% Get the curve fit result names
curve_fit_results_names = fieldnamesr(curve_fit_results,'prefix');

% Plot the ttest results
% plot_ttest_results(curve_fit_results, curve_fit_results_names, curve_fit_mode, ttest_results)
end

%% Curve fit function
function result = curve_fit_ikaro(result_matrix_3d, curve_fit_mode)

if strcmp(curve_fit_mode,'diff')    % If the 'fit' is going to be done be diff

    dif_teste = diff(result_matrix_3d,1,1);     % Get the difference between hours
    sum_dif_teste = nanmean(dif_teste);          % Sum the differences to get the direction of the change
    result = permute(sum_dif_teste,[3 2 1]);    % Store it into the result variable

else
    % Pre-allocate the results
    result = nan(size(result_matrix_3d,3),size(result_matrix_3d,2));
    for ii = 1:size(result_matrix_3d,2)        % Each animal
        for jj = 1:size(result_matrix_3d,3)    % Each day
            % Get the values regarding this specific parameter
            y_vector = result_matrix_3d(:,ii,jj);
            % Get the x vector (number of hours! important, since different
            % animals/days might have blank values
            x_vector = find(~isnan(y_vector));

            % Exclude the 7th hour
            x_vector(x_vector == 7) = [];

            if strcmp(curve_fit_mode,'fit')
                coefficients = polyfit(x_vector, y_vector(x_vector), 1);
                % Extract the slope (linear coefficient);
                linear_coefficient = coefficients(1);
                % Get the results
                result(jj,ii) = linear_coefficient;
            elseif strcmp(curve_fit_mode,'corr')
                R = corr(x_vector, y_vector(x_vector));
                result(jj,ii) = R;
            end
        end
    end

end
end

%% Calculate anova

function [anova_results,anova_rejected] = calculate_anova(results)
% Get every single results field
results_names = fieldnamesr(results,'prefix');
% Create a table for anova results
anova_results = cell(length(results_names),1);

for ii = 1:length(results_names)

    % execute command
    eval(sprintf("[anova_results{ii,1},anova_results{ii,2},anova_results{ii,3}] = anova1(%s',[],'off');",results_names{ii}))

end

% Get the results names
anova_results = [anova_results results_names];

% Get the comparisons which wielded a p-value < 0.05
rejected = find([anova_results{:,1}] < 0.05);
anova_rejected = {results_names{rejected}; anova_results{rejected,1}; anova_results{rejected,2}; anova_results{rejected,3}}';
end

%% Plot function
function plot_curve_fit_multiple_days(results, result_name, curve_fit_mode)
default_colors = lines;
save_coeff = [];
figure
for ii = 1:size(results,2)        % Each animal
    hold on
    for jj = 1:size(results,3)    % Each day
        % Get the values regarding this specific parameter
        y_vector = results(:,ii,jj);
        % Get the x vector (number of hours! important, since different
        % animals/days might have blank values
        x_vector = find(~isnan(y_vector));

        % Exclude the 7th hour
        x_vector(x_vector == 7) = [];

        coefficients = polyfit(x_vector, y_vector(x_vector), 1);

        add_jump = 8*(jj-1);

        subplot(3,1,1)
        title(result_name)
        % Plot the scatter
        scatter(add_jump+x_vector, y_vector(x_vector), 'o', 'filled','MarkerFaceColor',default_colors(ii,:));  % Scatter plot of the original data points
        hold on;

        % Plot the linear fit line
        fit_line = polyval(coefficients, x_vector);
        plot(add_jump+x_vector, fit_line, 'LineWidth', 2, 'Color',default_colors(ii,:));  % Plot the linear fit line in red
        xticks([])


        subplot(3,1,2)
        coefficients(2) = 0;

        if strcmp(curve_fit_mode,'fit')
            save_coeff(jj,ii) = coefficients(1);
        elseif strcmp(curve_fit_mode,'corr')
            R = corr(x_vector, y_vector(x_vector));
            save_coeff(jj,ii) = R;
            coefficients(1) = R;
        end

        % Fit the line
        fit_line = polyval(coefficients, x_vector);
        plot(add_jump+x_vector, fit_line, 'LineWidth', 2, 'Color',default_colors(ii,:));  % Plot the linear fit line in red
        xticklabels({})

    end
end

% Plot the average linear coeff as a black line
av_coeff = mean(save_coeff');
sem_coeff = std(save_coeff')/sqrt(size(save_coeff,2));
x_vector = 0:7;
subplot(3,1,2)
hold on
for ii = 1:size(results,3)  % Plot an average for each day
    add_jump = 8*(ii-1);
    % Fit the line
    fit_line = polyval([av_coeff(ii) 0], x_vector);
    plot(add_jump+x_vector, fit_line, 'LineWidth', 3, 'Color','k');  % Plot the linear fit line in red
    xticks([])
    yticks([])

    % Add the angle at the end of the black line
    %     text(add_jump+x_vector(end),fit_line(end),string(atan(av_coeff(ii))))

end

% Plot a bar plot showing the average
subplot(3,1,3)
bar(categorical({'BL','D1','D2','D3','D4'}),av_coeff)    % Plot averages
if strcmp(curve_fit_mode,'fit')
    ylabel('Av. angular coeff')
elseif strcmp(curve_fit_mode,'corr')
    ylabel('Corr Coeff')
end

hold on
for ii = 1:5
    plot([ii ii],[av_coeff(ii)-sem_coeff(ii) av_coeff(ii)+sem_coeff(ii)],'Color','k','LineWidth',2)
    plot(ii,save_coeff(ii,:),'o','Color','k')
end

% saveas(gcf,fullfile('E:\Barnes Maze - Mestrad\dados matlab\CurveFitResults',sprintf('%s.eps',result_name)))
% Save the figures
if strcmp(curve_fit_mode,'fit')
    saveas(gcf,fullfile('E:\Barnes Maze - Mestrad\dados matlab\CurveFitResults',sprintf('%s.png',result_name)))
    pause(0.2)
elseif strcmp(curve_fit_mode,'corr')
    saveas(gcf,fullfile('E:\Barnes Maze - Mestrad\dados matlab\CurveCorrResults',sprintf('%s.png',result_name)))
    pause(0.2)
end

close gcf
end

%% Calculate ttest2

function [ttest_results, filtered_ttest] = run_ttest(results)
% Get every single results field
results_names = fieldnamesr(results,'prefix');
% Create a table for anova results
ttest2_results = cell(length(results_names),1);

for ii = 1:length(results_names)
    eval(sprintf("matrix = %s';",results_names{ii}));
    for jj = 1:size(matrix,2)
        for kk = 1:size(matrix,2)
            eval(sprintf("a = %s';",results_names{ii}))
            a = a(:,jj);
            eval(sprintf("b = %s';",results_names{ii}))
            b = b(:,kk);
            % execute command
            [results_tt_h(jj,kk), results_tt_p(jj,kk)] = ttest(a,b);
            %eval(sprintf("[ttest2_results{ii,1},ttest2_results{ii,2}] = ttest2(%s',[],'off');",results_names{ii}))
        end
    end
    ttest_results{ii,1} = results_tt_h;
    ttest_results{ii,2} = results_tt_p;
    ttest_results{ii,3} = results_names{ii};
    clear results_tt_h results_tt_p
end

% Filter results
filtered_ttest = {};
count = 1;
for ii = 1:size(ttest_results,1)
    if nansum(find(ttest_results{ii,1} == 1)) > 1
        filtered_ttest{count,1} = ttest_results{ii,1};
        filtered_ttest{count,2} = ttest_results{ii,2};
        filtered_ttest{count,3} = ttest_results{ii,3};
        count = count + 1;
    end
end

end


%% Plot tttest results

function plot_ttest_results(curve_fit_results, results_names, curve_fit_mode, ttest_results)
default_colors = lines;
save_coeff = [];
figure
for ii = 1:size(ttest_results,1)
    for jj = 1:size(curve_fit_results,1)    % Each result
        t_matrix = ttest_results{jj};
        r_name = results_names{jj};
        eval(sprintf('result = %s;',results_names{jj}));

        % Exclude the 7th hour
        x_vector = 1:6;

        for kk = 1:size(result,2)
            coefficients = [result(jj,kk) 0];

            add_jump = 8*(jj-1);

            subplot(2,1,1)
            title(r_name)

            % Plot the linear fit line
            fit_line = polyval(coefficients, x_vector);
            hold on
            plot(add_jump+x_vector, fit_line, 'LineWidth', 2, 'Color',default_colors(kk,:));  % Plot the linear fit line in red
            xticks([])

        end


        % Plot the average linear coeff as a black line
        av_coeff = mean(result(jj,:));
        sem_coeff = std(result(jj,:))/sqrt(length(result(jj,:)));
        x_vector = 0:7;
        subplot(2,1,1)
        hold on

        % Fit the line
        fit_line = polyval([av_coeff 0], x_vector);
        plot(add_jump+x_vector, fit_line, 'LineWidth', 3, 'Color','k');  % Plot the linear fit line in red
        xticks([])
        yticks([])


    end

    % Plot a bar plot showing the average
    subplot(2,1,2)
    av_coeff = mean(result');
    sem_coeff = std(result')/sqrt(length(result(1,:)));
    bar(categorical({'BL','D1','D2','D3','D4'}),av_coeff)    % Plot averages
    if strcmp(curve_fit_mode,'fit')
        ylabel('Av. angular coeff')
    elseif strcmp(curve_fit_mode,'corr')
        ylabel('Corr Coeff')
    end

    hold on
    for kk = 1:5
        plot([kk kk],[av_coeff(kk)-sem_coeff(kk) av_coeff(kk)+sem_coeff(kk)],'Color','k','LineWidth',2)
        plot(kk,result(kk,:),'o','Color','k')
    end

    % Plot the ttest results
end



if strcmp(curve_fit_mode,'fit')
%     saveas(gcf,fullfile('E:\Barnes Maze - Mestrad\dados matlab\CurveFitResults',sprintf('%s.png',results_names)))
    saveas(gcf,fullfile('E:\Barnes Maze - Mestrad\dados matlab\CurveFitResults',sprintf('%s.eps',result_name)))
    % Save the figures
    pause(0.2)
elseif strcmp(curve_fit_mode,'corr')
    saveas(gcf,fullfile('E:\Barnes Maze - Mestrad\dados matlab\CurveCorrResults',sprintf('%s.png',results_names)))
    pause(0.2)
end

close gcf
end