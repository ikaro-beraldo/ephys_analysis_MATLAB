function create_prism_files_for_results_first_day(results_names, normalized_results, dir_to_save)


model_filename = "E:\Barnes Maze - Mestrad\dados matlab\PrismResults\ModelFile - First Day.pzfx";


% Lop for each result
for r = 1:length(results_names)
    % Define the command
    command = sprintf("new_matrix = normalized_%s;",results_names{r});
    % Execute the command
    eval(command)
    % Permute the results
    new_matrix = permute(new_matrix,[3 2 1]);
    % Take the Time_bin results and extract the first day (transpose it so
    % the rows = subjects and columns = days (prism default)
    new_matrix = new_matrix(:,:,1)';
    % Get the result name (also change the dots (.) to underscore (_)
    result_name = strrep(results_names{r},'.','_');
    % Define a final filename
    final_file = fullfile(dir_to_save, result_name);

    % Finally, create the pzfx (prism file)
    create_new_pzfx_file(model_filename,new_matrix,result_name,final_file)
    fprintf('%d - %s: OK!\n',r,results_names{r})
end
end