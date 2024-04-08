% Create prism files for results using two factor analysis (i.e. two way
% anova)

function create_prism_files_for_KW_analysis(results_names, results_distribution_bin, dir_to_save,flag_days, limit_hour)

% Analisys with 5 days
if flag_days == 5
    model_filename = "E:\Barnes Maze - Mestrad\dados matlab\PrismResultsKW\MoldeKW.pzfx";
elseif flag_days == 4   % Analysis with only 4 days
    model_filename = "E:\Barnes Maze - Mestrad\dados matlab\PrismResults\ModelFile2.pzfx";
elseif flag_days == 0   % non-parametric test
    model_filename = "E:\Barnes Maze - Mestrad\dados matlab\PrismResults\Friedman_molde.pzfx";
end

% Lop for each result
for r = 1:length(results_names)
    % Define the command
    command = sprintf("new_cell = %s;",results_names{r});
    % Execute the command
    eval(command)

    % IMPORTANT 2 --> transpose a matrix para hora virar dia e dia virar
    % hora
    new_cell = new_cell';

    % IMPORTANT 3 DAQUENAIPE PQ É BO QUE APARECE POR AQUI -->
    % LIMIT THE NUMBER OF HOURS BY EXCLUDING THE EXTRA ROWS
    if size(new_cell,1) > limit_hour
        new_cell(limit_hour+1:end,:) = [];
    end

    % Get the result name (also change the dots (.) to underscore (_)
    result_name = strrep(results_names{r},'.','_');
    % Define a final filename
    final_file = fullfile(dir_to_save, result_name);

    % Finally, create the pzfx (prism file)
    create_new_pzfx_file(model_filename,new_cell,result_name,final_file)
    fprintf('%d - %s: OK!\n',r,results_names{r})
end
end

%% Helper functions

%% Create a new pzfx file (Prism file) based on a model

function create_new_pzfx_file(model_filename,new_cell,result_name,final_file)
% Load the prism model
model = importdata(model_filename);

%Define which line to be changed
rows_to_change = find(strcmp(model,'<d>0</d>'));

% Get the total number of elements that will be changed
n_elements = 0;
for ii = 1:numel(new_cell)
    n_elements = n_elements + length(new_cell{ii});
end

% Subtract the ones that will not be added 
n_elements = n_elements - length(rows_to_change);

% Create a final model (considering the empty elements on the cell that
% will be filled
final_model = cell(length(model) + n_elements,1);

% Add the model elements to the final_model
final_model(1:length(model),1) = model(:,1);

% Now for the changes!
% Change the Parameter name (line 23)
final_model{24,1} = sprintf('<Title>%s</Title>',result_name);

% Counter to be add to the cell index
added_idx = 0;
num_elements_to_add = 1;

for i = 1:length(rows_to_change)

    % Increase the index after adding new elements
    rows_to_change = rows_to_change + num_elements_to_add-1;
    % Get the current index
    index = rows_to_change(i);
      
    % Add elements after the current index
    num_elements_to_add = length(new_cell{i});

    % Find the last filled cell (to get which is the last element)
    filled_cells = find(~cellfun(@isempty,final_model));
    save_elements = final_model(index+1:filled_cells(end));

    % Save the elements
    final_model((index+1:filled_cells(end))+num_elements_to_add-1) = save_elements;

    % Change the selected index itself (before adding new elements)
    final_model{index} = sprintf('<d>%d</d>',new_cell{i}(1));

    % Loop for each single value added
    for j = 2:num_elements_to_add
        % Get the value that must be added
        added_value = sprintf('<d>%d</d>',new_cell{i}(j));

        % Add the value
        final_model{index+1} = added_value;       
        
        index = index + 1;
    end
end


% Create a txt file with the new info
txt_final_file = sprintf('%s.txt',final_file);
% Create a txt file with the new infotxt_final_file = sprintf('%s.txt',final_file);
writecell(final_model,txt_final_file,'QuoteStrings','none')
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