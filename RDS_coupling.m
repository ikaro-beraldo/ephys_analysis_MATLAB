%% Function to the get the Ripple-Delta-Spindle coupling
function [linear, block, number] = RDS_coupling(ripple_blocks, delta_blocks, spindle_blocks, fs, block_length)

%% Event checker
ev_checker = [false false false]; %[ripple spindle delta]

%% Transform the block indexing into vector indexing (for comparison)
spindle_timestamp = fix_block_to_linear(spindle_blocks,fs, block_length);
delta_timestamp = fix_block_to_linear(delta_blocks,fs, block_length);
ripple_timestamp = fix_block_to_linear(ripple_blocks,fs, block_length);

if isempty(ripple_timestamp) == false
    ev_checker(1) = true;
end
if isempty(spindle_timestamp) == false
    ev_checker(2) = true;
end
if isempty(delta_timestamp) == false
    ev_checker(3) = true;
end

%% Get event coupling

%Delta - spindle
if ev_checker(3) == true && ev_checker(2) == true
    peak_dist_min = 50 * fs/1000; %100
    peak_dist_max = 1500 * fs/1000; %1300

    delta_spindle = spindle_timestamp(:,2) - delta_timestamp(:,2)';
    [l,c] = find(delta_spindle >= peak_dist_min & delta_spindle <= peak_dist_max);
    clear delta_spindle
    delta_spindle = [l c];
else
    delta_spindle = [];
end

% Ripple - Delta
if ev_checker(1) == true && ev_checker(3) == true
    peak_dist_min = 50 * fs/1000; %50
    peak_dist_max = 300 * fs/1000; %250

    ripple_delta = delta_timestamp(:,2) - ripple_timestamp(:,2)';
    [l,c] = find(ripple_delta >= peak_dist_min & ripple_delta <= peak_dist_max);
    clear ripple_delta
    ripple_delta = [l c];

    %Delta - Ripple
    peak_dist_min = 50 * fs/1000; %50
    peak_dist_max = 600 * fs/1000; %400

    delta_ripple = ripple_timestamp(:,2) - delta_timestamp(:,2)';
    [l,c] = find(delta_ripple>= peak_dist_min & delta_ripple <= peak_dist_max);
    clear delta_ripple
    delta_ripple = [l c];
else
    ripple_delta = [];
    delta_ripple = [];
end

%Ripple - Delta - Spindle
if isempty(ripple_delta) == false && isempty(delta_spindle) == false
    ripple_delta_spindle = [];
    for rds = 1:size(ripple_delta,1)
        compare = ripple_delta(rds,1) == delta_spindle(:,2);
        if sum(compare) >= 1
            compare_idx = find(compare == 1);
            ripple_delta_spindle = [ripple_delta_spindle;...
                ripple_delta(rds,:) delta_spindle(compare_idx(1),1)];
        end

    end
else
    ripple_delta_spindle = [];
end
%% Number of events
number.delta_ripple = size(delta_ripple,1);
number.ripple_delta = size(ripple_delta,1);
number.delta_spindle = size(delta_spindle,1);
number.ripple_delta_spindle = size(ripple_delta_spindle,1);

%% Timestamps

if isempty(delta_ripple) == false
    linear.delta_ripple_timestamps.delta = delta_timestamp(delta_ripple(:,2),:);
    linear.delta_ripple_timestamps.ripple = ripple_timestamp(delta_ripple(:,1),1:3);
end

if isempty(ripple_delta) == false
    linear.ripple_delta_timestamps.ripple = ripple_timestamp(ripple_delta(:,2),:);
    linear.ripple_delta_timestamps.delta = delta_timestamp(ripple_delta(:,1),:);
end

if isempty(delta_spindle) == false
    linear.delta_spindle_timestamps.delta = delta_timestamp(delta_spindle(:,2),:);
    linear.delta_spindle_timestamps.spindle = spindle_timestamp(delta_spindle(:,1),:);
end

if isempty(ripple_delta_spindle) == false
    linear.ripple_delta_spindle_timestamps.ripple = ripple_timestamp(ripple_delta_spindle(:,2),:);
    linear.ripple_delta_spindle_timestamps.delta = delta_timestamp(ripple_delta_spindle(:,1),:);
    linear.ripple_delta_spindle_timestamps.spindle = spindle_timestamp(ripple_delta_spindle(:,3),:);
end

%% Get the block indices

if isempty(delta_ripple) == false
    block.delta_ripple_timestamps.delta = fix_linear_to_block(linear.delta_ripple_timestamps.delta,fs, block_length);
    block.delta_ripple_timestamps.ripple = fix_linear_to_block(linear.delta_ripple_timestamps.ripple,fs, block_length);
end

if isempty(ripple_delta) == false
    block.ripple_delta_timestamps.ripple = fix_linear_to_block(linear.ripple_delta_timestamps.ripple,fs, block_length);
    block.ripple_delta_timestamps.delta = fix_linear_to_block(linear.ripple_delta_timestamps.delta,fs, block_length);
end

if isempty(delta_spindle) == false
    block.delta_spindle_timestamps.delta = fix_linear_to_block(linear.delta_spindle_timestamps.delta,fs, block_length);
    block.delta_spindle_timestamps.spindle = fix_linear_to_block(linear.delta_spindle_timestamps.spindle,fs, block_length);
end

if isempty(ripple_delta_spindle) == false
    block.ripple_delta_spindle_timestamps.ripple = fix_linear_to_block(linear.ripple_delta_spindle_timestamps.ripple,fs, block_length);
    block.ripple_delta_spindle_timestamps.delta = fix_linear_to_block(linear.ripple_delta_spindle_timestamps.delta,fs, block_length);
    block.ripple_delta_spindle_timestamps.spindle = fix_linear_to_block(linear.ripple_delta_spindle_timestamps.spindle,fs, block_length);
end

end

%% Helper function
function linear_result = fix_block_to_linear(result_block,fs, block_length)
block_length_total = fs * block_length;
factor = block_length_total * (result_block(:,1)-1);

% Add the factor to the rest
linear_result = result_block(:,2:4)+factor;
end

function block_result = fix_linear_to_block(result_linear,fs,block_length)
block_length_total = fs * block_length;

block_result = [ceil(result_linear(:,1)/block_length_total) mod(result_linear,block_length_total)];
end