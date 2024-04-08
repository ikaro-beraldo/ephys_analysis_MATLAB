%% Create a new pzfx file (Prism file) based on a model

function create_new_pzfx_file(model_filename,new_matrix,result_name,final_file)
% Load the prism model
model = importdata(model_filename);

%Define which line to be changed
rows_to_change = find(strcmp(model,'<d>1</d>'));

% Now for the changes!
% Change the Paramete name (line 23)
model{23,1} = sprintf('<Title>%s</Title>',result_name);

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

%%


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