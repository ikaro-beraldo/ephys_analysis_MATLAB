clear all

frame_rate = 20;    % Frame rate não esquece mateus
% CSV file name
fileName = 'C:\Users\ikaro\Downloads\teste2.csv'; % CSV path\filename
% Open CSV data as a table
T = readtable(fileName);
% Get the Behaviour names
vid_num = T.x_CSV_HEADER_Metadata_id;
behav_names = T.metadata;
behav_start = T.temporal_segment_start; % Start
behav_end = T.temporal_segment_end;     % End

% for loop 
for ii = 1:length(behav_names)
    a = behav_names{ii};    % selected one
    b = a(23:end-2);        % Filter unimportant information
    behav_names{ii} = b;

    % Filt the video number
    vid_num{ii} = str2double(vid_num{ii}(1));
end

% Get the video_numbers
vid_num = cell2mat(vid_num);
% Get the unique behavior names
behav_list = unique(behav_names);

%% Separate videos one by one
% Get the unique video numbers
all_vid_num = unique(vid_num);
behav_cell_final = cell(length(all_vid_num),2); % Cell with final results

% Loop for each video
for ii = 1:length(all_vid_num)
    vid_idx = all_vid_num(ii);
    % Extract only the rows regarding the selected video
    selected_vid = find(vid_num == vid_idx);
    s_behav_end = behav_end(selected_vid);
    s_behav_start = behav_start(selected_vid);
    s_behav_names = behav_names(selected_vid);

    % Get the unique behaviours on this specific video
    unique_vid_behav_names = unique(s_behav_names);
    behav_cell_final{ii,2} = unique_vid_behav_names;

    time_each_behaviour = nan(length(unique_vid_behav_names),3);
    % Logical matrix with
    behav_matrix = zeros(length(unique_vid_behav_names),round(s_behav_end(end)*frame_rate)+1);

    % Loop for each behaviour
    for jj = 1:length(unique_vid_behav_names)
        selected_behav = unique_vid_behav_names{jj};    % selected behaviour
        % Get the moments at which the behavior is true
        moments = strcmp(selected_behav,s_behav_names);
        
        % Convert to sec2frame
        starts = round(s_behav_start(moments)*frame_rate)+1;
        finals = round(s_behav_end(moments)*frame_rate)+1;

        % Insert the values on the final matrix
        for kk = 1:length(starts)
            behav_matrix(jj,starts(kk):finals(kk)) = 1;
        end
        % Calculate time for each behaviour
        time_each_behaviour(jj,1) = sum(behav_matrix(jj,:));                        % In Frame
        time_each_behaviour(jj,2) = sum(behav_matrix(jj,:))/frame_rate;             % In seconds
    end
    
    time_each_behaviour(:,3) = sum(behav_matrix,2)/length(find(behav_matrix == 1));   % Percentage

    % Insert data on the FINAL cell 
    behav_cell_final{ii,1} = behav_matrix;
    behav_cell_final{ii,3} = time_each_behaviour;
end

% Output file
output_file = sprintf('%s.mat',fileName(1:end-4));
save(output_file,'behav_cell_final','frame_rate','T')

