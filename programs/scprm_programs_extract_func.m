% Extract veridical direction (a_unique), stimulus direction (a_stim) and
% coherence value (coh_val)

function [a_unique,a_stim,coh_val] = scprm_program_extract_func(bl_unique,dc,dc_mat,dc_ts,rdp_dir,d_rdp_dir_ts,dp_all,dp_all_ts)

coh_val = dc{bl_unique};
for i = 1:length(bl_unique)
    if bl_unique(i) == length(dc_mat)
        bl = d_rdp_dir_ts>=dc_ts(bl_unique(i));
    else
        bl = and(d_rdp_dir_ts>=dc_ts(bl_unique(i)),d_rdp_dir_ts< dc_ts(bl_unique(i)+1));
    end
    ss{i} = find(bl == 1); % indices steady states in that block
    nstate_cell{i} = length(ss{i});
    rdp_ss_cell{i} = deg2rad(cell2mat(rdp_dir(ss{i})));
    if bl_unique(i) == length(dc_mat)
        substate = and(dp_all_ts>=d_rdp_dir_ts(ss{i}(1)),dp_all_ts<d_rdp_dir_ts(ss{i}(end)));
    else
        substate = and(dp_all_ts>=d_rdp_dir_ts(ss{i}(1)),dp_all_ts<d_rdp_dir_ts(ss{i}(end)+1));
    end
    %substate = and(dp_all_ts>=d_rdp_dir_ts(ss{i}(1)),dp_all_ts<d_rdp_dir_ts(ss{i}(end)+1));
    ss_ts_cell{i} = find(substate == 1);
end
nstate = sum(cell2mat(nstate_cell));
n_life = 25;
ss_ts = cell2mat(ss_ts_cell);
dp_ts = dp_all_ts(ss_ts);
dp = dp_all(ss_ts);
rdp_ss_mat = cell2mat(rdp_ss_cell);

t_index = cell(length(nstate_cell),10); % maximum steady state is 10
for j = 1:length(nstate_cell)
    for i = 1:nstate_cell{j}
        if ss{j}(i) == length(rdp_dir)
            st = and(dp_all_ts>=d_rdp_dir_ts(ss{j}(i)),dp_all_ts<d_rdp_dir_ts(ss{j}(end)));
        else
            st = and(dp_all_ts>=d_rdp_dir_ts(ss{j}(i)),dp_all_ts<d_rdp_dir_ts(ss{j}(i)+1));
        end
        t_index{j,i} = find(st == 1); 
    end
end

a_unique = [];
for i = 1:length(nstate_cell)
    for j = 1:length(rdp_ss_cell{i})
        a_unique = [a_unique rdp_ss_cell{i}(j)*ones(1,length(t_index{i,j}))];
    end
end

a_unique(end) = [];

f = 1:length(dp);
for iFrame = 1:length(f)
    xIdx            = logical(mod([1:size(dp{iFrame},2)],2));
    xpos{iFrame}    = dp{iFrame}(xIdx);
    ypos{iFrame}    = dp{iFrame}(~xIdx);
end

ndots = length(xpos{1});

for i = 1:ndots
    for j = 1:length(f)-1
        x_resultant{i}(j) = xpos{j+1}(i) - xpos{j}(i);
        y_resultant{i}(j) = ypos{j+1}(i) - ypos{j}(i);
        dir{i}(j) = mod(atan2(x_resultant{i}(j),y_resultant{i}(j)),2*pi);
        d_length(i,j) = sqrt((x_resultant{i}(j))^2 + (y_resultant{i}(j))^2);
    end
end

a_true = zeros(length(dir),length(dir{1}));

for i = 1:length(dir) % direction
    for j = 1:length(dir{1}) % frame
        a_true(i,j) = dir{i}(j);
    end
end

dir_test = dir;


% Remove dots that disappear from d_length
d_length_mean = mean(mode(d_length));
for i = 1:ndots
    for j = 1:iFrame-1
        if abs(d_length(i,j)-d_length_mean) > 4.9e-4
            dir{i}(j) = NaN;
            x_resultant{i}(j) = NaN;
            y_resultant{i}(j) = NaN;
            a_true(i,j) = NaN;
        end
    end
end

for i = 1:ndots
    x_speed(i,:) = x_resultant{i};
    y_speed(i,:) = y_resultant{i};
end

for j = 1:iFrame-1
    ndots_avail(j) = ndots - length(find(isnan(x_speed(:,j))==1));
end
    

for i = 1:iFrame-1
    a_sig{i} = find(abs(circ_dist(a_true(:,i),a_unique(i))) < 1e-3);
    a_noise{i} = find(abs(circ_dist(a_true(:,i),a_unique(i))) >= 1e-3);
    x_sig{i} = x_speed(a_sig{i},i);
    x_noise{i} = x_speed(a_noise{i},i);
    y_sig{i} = y_speed(a_sig{i},i);
    y_noise{i} = y_speed(a_noise{i},i);
    coh_est(i) = length(a_sig{i})/ndots_avail(i);
end

for i = 1:iFrame-1
    x_sig_sum(i) = sum(x_sig{i}); 
    y_sig_sum(i) = sum(y_sig{i});
    x_noise_sum(i) = sum(x_noise{i});
    y_noise_sum(i) = sum(y_noise{i});
end

x_sig_mean = x_sig_sum./ndots_avail;
y_sig_mean = y_sig_sum./ndots_avail;
x_noise_mean = x_noise_sum./ndots_avail;
y_noise_mean = y_noise_sum./ndots_avail;

a_stim = mod(atan2(x_sig_mean+x_noise_mean,y_sig_mean+y_noise_mean),2*pi);
