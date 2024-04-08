% Create prism files for results using two factor analysis (i.e. two way
% anova)

function create_prism_files_for_two_factor_analysis(results_names, normalized_results, dir_to_save,flag_days, limit_hour)

% Analisys with 5 days
if flag_days == 5
    model_filename = "E:\Barnes Maze - Mestrad\dados matlab\PrismResults_hour-hour\Molde_hour-hour.pzfx";
elseif flag_days == 4   % Analysis with only 4 days
    model_filename = "E:\Barnes Maze - Mestrad\dados matlab\PrismResults\ModelFile2.pzfx";
elseif flag_days == 0   % non-parametric test
    model_filename = "E:\Barnes Maze - Mestrad\dados matlab\PrismResults\Friedman_molde.pzfx";
end

% Lop for each result
for r = 1:length(results_names)
    % Define the command
    command = sprintf("new_matrix = normalized_%s;",results_names{r});
    % Execute the command
    eval(command)

    % IMPORTANT --> CHECK IF THE RESULT IS A 3D ARRAY (CONSIDERING THE 2
    % FACTOR ANALYSIS)
    if ndims(new_matrix) ~= 3
        continue
    end

    % IMPORTANT 2 KKKKKK --> PERMUTE THE MATRIX, SO IT CAN BE EASILY
    % INSERTED INTO PRISM
    new_matrix = permute(new_matrix,[3 2 1]);

    % IMPORTANT 3 DAQUENAIPE PQ É BO QUE APARECE POR AQUI -->
    % LIMIT THE NUMBER OF HOURS BY EXCLUDING THE EXTRA ROWS
    if size(new_matrix,1) > limit_hour
        new_matrix(limit_hour+1:end,:,:) = [];
    end

    % Get the result name (also change the dots (.) to underscore (_)
    result_name = strrep(results_names{r},'.','_');
    % Define a final filename
    final_file = fullfile(dir_to_save, result_name);

    % Finally, create the pzfx (prism file)
    create_new_pzfx_file(model_filename,new_matrix,result_name,final_file)
    fprintf('%d - %s: OK!\n',r,results_names{r})
end
end

%% Helper functions

%% Create a new pzfx file (Prism file) based on a model

function create_new_pzfx_file(model_filename,new_matrix,result_name,final_file)
% Load the prism model
model = importdata(model_filename);

%Define which line to be changed
rows_to_change = find(strcmp(model,'<d>0</d>'));

% Now for the changes!
% Change the Parameter name (line 23)
model{24,1} = sprintf('<Title>%s</Title>',result_name);

% Change the values according to the input matrix
for ii = 1:length(rows_to_change)
    model{rows_to_change(ii),1} = sprintf('<d>%d</d>',new_matrix(ii));
end

% Create a txt file with the new info
txt_final_file = sprintf('%s.txt',final_file);
writecell(model,txt_final_file,'QuoteStrings','none')
% Create a copy with the extension pzfx
change_file_extension_to_pzfx(txt_final_file)
end

%% Function to change the file extension

function change_file_extension_to_pzfx(filename)
% Create sample filename.
[folder, baseFileName, ~] = fileparts(filename);
% Ignore extension and replace it with .txt
newBaseFileName = sprintf('%s.pzfx', baseFileName);
% Make sure the folder is prepended (if it has one).
newFullFileName = fullfile(folder, newBaseFileName);

% Now copy the content
copyfile(filename,newFullFileName)
end