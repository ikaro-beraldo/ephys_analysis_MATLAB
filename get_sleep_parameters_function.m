% Get sleep parameters for each day/animal

function get_sleep_parameters_function(directory,names)

s_p = Signal_Processing;

mkdir(fullfile(directory,'sleep_parameters'))
for ii = 1:length(names)

    clearvars -except directory names counter ii s_p
    counter = 1;

    % Complete filenames
    filename2 = fullfile(directory,names{ii},'GMM_Classification.mat');
    load(filename2)
    % Fix All_Sort
    GMM_NREM_All_Sort(:) = false;
    GMM_NREM_All_Sort(GMM.All_Sort == 2) = true;
    GMM_REM_All_Sort(:) = false;
    GMM_REM_All_Sort(GMM.All_Sort == 1) = true;
    GMM_WK_All_Sort(:) = false;
    GMM_WK_All_Sort(GMM.All_Sort == 3) = true;

    % Call function
    sleep_param = sleep_parameters(GMM.All_Sort);

    sav = fullfile(directory,'sleep_parameters',char(names{ii}));
    save(sav,'sleep_param')

    disp(names{ii})

end
end