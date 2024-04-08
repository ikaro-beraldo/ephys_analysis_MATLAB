% Detect RDS block by block
function detect_RDS_block_by_block_function(directory,names,names2)

mkdir(fullfile(directory,'RDS'))
for ii = 1:length(names)
    disp(names{ii})
    clearvars -except directory names* counter ii s_p
    counter = 1;
    
    % Complete filename
    filename1 = fullfile(directory,'ripple',names{ii});
    filename2 = fullfile(directory,'delta',names{ii});
    filename3 = fullfile(directory,'spindle',names{ii});
    filename4 = fullfile(directory,names{ii},'GMM_Classification.mat');
    filename5 = fullfile(directory,names{ii},'blocked_data.mat');
    load(filename1,'ripple_result_block_2','ripple_result_block_3','ripple_accepted')
    load(filename2,'delta_blocks')
    load(filename3,'spindle_blocks')
    load(filename5,'LFP1','fs')
    load(filename4,'GMM')
    
    % Fix Delta blocks (since it didn't considered all the epochs (only
    % NREM)
    % nrem = find(GMM.All_Sort == 2);
    % delta_blocks(:,1) = nrem(delta_blocks(:,1));
    block_length = size(LFP1,2)/fs;
    
    % Detect RDS coupling
    [linear_2, block_2, number_2] = RDS_coupling(ripple_result_block_2, delta_blocks, spindle_blocks, fs, block_length);
    % Detect RDS coupling
    [linear_3, block_3, number_3] = RDS_coupling(ripple_result_block_3, delta_blocks, spindle_blocks, fs, block_length);
    % Detect RDS coupling
    [linear, block, number] = RDS_coupling(ripple_accepted, delta_blocks, spindle_blocks, fs, block_length);

    
    %display(['total number of candidate events are ' num2str(size(spindles_blocks,1))])
    
    sav = fullfile(directory,'RDS',names{ii});
    save(sav,'linear','block','number','linear_2','block_2','number_2','linear_3','block_3','number_3')
    
    
end

%% Get the number of events and coupling for each animal/day
block_length = 10;

for ii = 1:size(names2,1)
    for jj = 1:size(names2,2)
        disp(names2{ii,jj})
        clearvars -except directory names* counter ii jj s_p ripple* delta* spindle* block_length absolute relative
        counter = 1;
        
        % Complete filename
        filename1 = fullfile(directory,'ripple',names2{ii,jj});
        filename2 = fullfile(directory,'delta',names2{ii,jj});
        filename3 = fullfile(directory,'spindle',names2{ii,jj});
        filename4 = fullfile(directory,'RDS',names2{ii,jj});
        filename5 = fullfile(directory,names2{ii,jj},'GMM_Classification.mat');
        
        load(filename1,'ripple_result_block_2','ripple_result_block_3','ripple_accepted')
        load(filename2,'delta_blocks')
        load(filename3,'spindle_blocks')
        load(filename4,'number_2','number_3','number')
        load(filename5,'GMM')
        
        % Get the number of isolated events
        absolute.ripple_2(ii,jj) = size(ripple_result_block_2,1);
        absolute.ripple_3(ii,jj) = size(ripple_result_block_3,1);
        absolute.ripple(ii,jj) = size(ripple_accepted,1);
        absolute.delta(ii,jj) = size(delta_blocks,1);
        absolute.spindle(ii,jj) = size(spindle_blocks,1);
        % get the number of coupled events Channel 2
        absolute.ripple_delta_2(ii,jj) = number_2.ripple_delta;
        absolute.delta_ripple_2(ii,jj) = number_2.delta_ripple;
        absolute.delta_spindle(ii,jj) = number_2.delta_spindle;
        absolute.ripple_delta_spindle_2(ii,jj) = number_2.ripple_delta_spindle;
        % get the number of coupled events Channel 3
        absolute.ripple_delta_3(ii,jj) = number_3.ripple_delta;
        absolute.delta_ripple_3(ii,jj) = number_3.delta_ripple;
        absolute.ripple_delta_spindle_3(ii,jj) = number_3.ripple_delta_spindle;
        % get the number of coupled events Accepted Ripples
        absolute.ripple_delta(ii,jj) = number.ripple_delta;
        absolute.delta_ripple(ii,jj) = number.delta_ripple;
        absolute.ripple_delta_spindle(ii,jj) = number.ripple_delta_spindle;
        
        nrem_length_minutes = length(find(GMM.All_Sort == 2))*block_length;
        % Get the number of isolated events
        relative.ripple_2(ii,jj) = size(ripple_result_block_2,1)/nrem_length_minutes;
        relative.ripple_3(ii,jj) = size(ripple_result_block_3,1)/nrem_length_minutes;
        relative.ripple(ii,jj) = size(ripple_accepted,1)/nrem_length_minutes;
        relative.delta(ii,jj) = size(delta_blocks,1)/nrem_length_minutes;
        relative.spindle(ii,jj) = size(spindle_blocks,1)/nrem_length_minutes;
        % get the number of coupled events Channel 2
        relative.ripple_delta_2(ii,jj) = number_2.ripple_delta/nrem_length_minutes;
        relative.delta_ripple_2(ii,jj) = number_2.delta_ripple/nrem_length_minutes;
        relative.delta_spindle(ii,jj) = number_2.delta_spindle/nrem_length_minutes;
        relative.ripple_delta_spindle_2(ii,jj) = number_2.ripple_delta_spindle/nrem_length_minutes;
        % get the number of coupled events Channel 2
        relative.ripple_delta_3(ii,jj) = number_3.ripple_delta/nrem_length_minutes;
        relative.delta_ripple_3(ii,jj) = number_3.delta_ripple/nrem_length_minutes;
        relative.ripple_delta_spindle_3(ii,jj) = number_3.ripple_delta_spindle/nrem_length_minutes;
        % get the number of coupled events Accepted Ripples
        relative.ripple_delta(ii,jj) = number.ripple_delta/nrem_length_minutes;
        relative.delta_ripple(ii,jj) = number.delta_ripple/nrem_length_minutes;
        relative.ripple_delta_spindle(ii,jj) = number.ripple_delta_spindle/nrem_length_minutes;
    end
end

% Save results
sav = fullfile(directory,'RDS','totalRDS');
save(sav,'relative','absolute')
end