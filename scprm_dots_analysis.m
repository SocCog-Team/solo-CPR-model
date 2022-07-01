% Analyze dots' data 
% IMPORTANT: x-axis is at 12 o'clock and y-axis is at 3 o'clock

clear all; close all;

addpath /Users/katenguyen/Dropbox/DPZ/MATLAB/Data

load('20220215_dem_CPRsolo_block2_fxs.mat') % load the mat file

idx.dot         = d.event == 'STIM_RDP_dot_positions';
dp_all          = d.value(idx.dot);
dp_all_ts       = d.time(idx.dot);

idx.coh         = d.event == 'STIM_RDP_coherence';
dc              = d.value(idx.coh);
dc_ts           = d.time(idx.coh);

idx.rdp_dir     = d.event == 'STIM_RDP_direction';
rdp_dir         = d.value(idx.rdp_dir);
d_rdp_dir_ts    = d.time(idx.rdp_dir); % frame time stamp

idx.jsDir       = d.event == 'IO_joystickDirection';
dp_dir          = d.value(idx.jsDir);
d_jsDir_ts      = d.time(idx.jsDir); % joystick time stamp

idx.jsStr       = d.event == 'IO_joystickStrength';
dp_str          = d.value(idx.jsStr);
d_jsStr_ts      = d.time(idx.jsStr); % joystick time stamp

bl = 4; % Choose your block
coh_val = dc{bl}; % coherence for that block
idx.bl = and(d_rdp_dir_ts>=dc_ts(bl),d_rdp_dir_ts< dc_ts(bl+1)); % indices of steady states within that block
ss = find(idx.bl == 1); % indices steady states in that block
nstate = length(ss);
substate = and(dp_all_ts>=d_rdp_dir_ts(ss(1)),dp_all_ts<d_rdp_dir_ts(ss(end)+1));
ss_ts = find(substate == 1);
dp_ts = dp_all_ts(ss_ts);
n_life = 25;

for i = 1:length(dp_ts)-1
    ti_diff(i) = 1e-6*(dp_ts(i+1) - dp_ts(i));
end

dp = dp_all(ss_ts);

f = 1:length(dp);
for iFrame = 1:length(f)
    xIdx            = logical(mod([1:size(dp{iFrame},2)],2));
    xpos{iFrame}    = dp{iFrame}(xIdx);
    ypos{iFrame}    = dp{iFrame}(~xIdx);
end

ndots = length(xpos{1});

for i = 1:ndots
    for j = 1:length(f)-1
        x_diff{i}(j) = xpos{j+1}(i) - xpos{j}(i); % x-distance that dot i travels between frame j and j+1
        y_diff{i}(j) = ypos{j+1}(i) - ypos{j}(i); % y-distance that dot i travels between frame j and j+1
        dir{i}(j) = mod(atan2(x_diff{i}(j),y_diff{i}(j)),2*pi); % direction that dot i travels between frame j and j+1
        d_length(i,j) = sqrt((x_diff{i}(j))^2 + (y_diff{i}(j))^2); % distance that dot i travels between frame j and j+1
    end
end

a_dot = zeros(length(dir),length(dir{1})); 

for i = 1:length(dir) % direction
    for j = 1:length(dir{1}) % frame
        a_dot(i,j) = dir{i}(j); % direction of dot i on frame j
    end
end


for i = 1:nstate
    st = and(dp_all_ts>=d_rdp_dir_ts(ss(i)),dp_all_ts<d_rdp_dir_ts(ss(i)+1));
    t_index{i} = find(st == 1); 
    a_veridical{i} = deg2rad(cell2mat(rdp_dir(ss(i))))*ones(1,length(t_index{i})); % veridical direction
end

a_unique_mat = cell2mat(a_unique); 

d_length_mean = mean(mode(d_length));
% Remove dots that has d_length > d_length_mean
for i = 1:ndots
    for j = 1:iFrame-1
        if abs(d_length(i,j)-d_length_mean) > 4.9e-4
            dir{i}(j) = NaN;
            x_resultant{i}(j) = NaN;
            y_resultant{i}(j) = NaN;
            a_dot(i,j) = NaN;
        end
    end
end

% Convert dots' distance between frames from cell to matrix
for i = 1:ndots
    x_dist(i,:) = x_diff{i}; % x-distance of dot i between frames
    y_dist(i,:) = y_diff{i}; % y-distance of dot i between frames
end

for j = 1:iFrame-1
    ndots_avail(j) = ndots - length(find(isnan(x_dist(:,j))==1)); % dots that don't disappear 
end

sum_x = sum(x_dist,'omitnan')./ndots_avail; % x-evidence per frame;
sum_y = sum(y_dist,'omitnan')./ndots_avail; % y-evidence per frame;


for i = 1:length(sum_x)
    dir_evidence(i) = mod(atan2(sum_x(i),sum_y(i)),2*pi);
end

for i = 1:iFrame-1
    a_sig{i} = find(abs(circ_dist(a_true(:,i),a_unique_mat(i))) < 1e-3);
    coh_est(i) = length(a_sig{i})/ndots_avail(i);
end

a_stimulus = mod(atan2(sum_x,sum_y),2*pi);