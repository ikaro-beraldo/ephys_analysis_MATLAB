function create_prism_files_for_results(results_names, normalized_results, dir_to_save,flag_days)

if flag_days == 5
    % Get prism file model
    model_filename = "E:\Barnes Maze - Mestrad\dados matlab\PrismResults\ModelFile - Reference.pzfx";
elseif flag_days == 4
    model_filename = "E:\Barnes Maze - Mestrad\dados matlab\PrismResults\ModelFile2.pzfx";
elseif flag_days == 0   % non-parametric test
    model_filename = "E:\Barnes Maze - Mestrad\dados matlab\PrismResults\Friedman_molde.pzfx";
end

% Lop for each result
for r = 1:length(results_names)
    % Define the command
    command = sprintf("new_matrix = normalized_%s';",results_names{r});
    % Execute the command
    eval(command)
    % Get the result name (also change the dots (.) to underscore (_)
    result_name = strrep(results_names{r},'.','_');
    % Define a final filename
    final_file = fullfile(dir_to_save, result_name);

    % Finally, create the pzfx (prism file)
    create_new_pzfx_file(model_filename,new_matrix,result_name,final_file)
    fprintf('%d - %s: OK!\n',r,results_names{r})
end
end