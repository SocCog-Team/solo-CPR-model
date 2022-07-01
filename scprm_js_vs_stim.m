% Compute optimal parameters to minimize the mean squared error between
% joystick direction and stimulus

clear all; close all;
%% At DPZ

addpath /Users/katenguyen/Dropbox/DPZ/MATLAB/Data
addpath /Users/katenguyen/Dropbox/DPZ/MATLAB/Programs

load('20220215_dem_CPRsolo_block2_fxs.mat') % load the mat file
%load('forKate.mat')

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
coh_val = dc{bl};
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


for i = 1:nstate
    st = and(dp_all_ts>=d_rdp_dir_ts(ss(i)),dp_all_ts<d_rdp_dir_ts(ss(i)+1));
    t_index{i} = find(st == 1); 
    a_unique{i} = deg2rad(cell2mat(rdp_dir(ss(i))))*ones(1,length(t_index{i}));
end

a_unique_mat = cell2mat(a_unique);

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

sum_x = sum(x_speed,'omitnan')./ndots_avail; % x-evidence per frame;
sum_y = sum(y_speed,'omitnan')./ndots_avail; % y-evidence per frame;


for i = 1:length(sum_x)
    dir_evidence(i) = mod(atan2(sum_x(i),sum_y(i)),2*pi);
end

for i = 1:iFrame-1
    a_sig{i} = find(abs(circ_dist(a_true(:,i),a_unique_mat(i))) < 1e-3);
    coh_est(i) = length(a_sig{i})/ndots_avail(i);
end

%% X_stim1
for i = 1:length(a_unique_mat)-1
    mu(:,i) = coh_val*0.0667*[sin(a_unique_mat(i));cos(a_unique_mat(i))];
end
%VX = 0.5*(0.0667)^2*(1-coh_val)/ndots;
VX = 0.0667^2*(1-mean(coh_est,'omitnan'))/(n_life*ndots);
S = [VX 0;0 VX];

for i = 1:length(mu)
    X_stim1(:,i) = mvnrnd(mu(:,i),S);
end

dir_stim1 = mod(atan2(X_stim1(1,:),X_stim1(2,:)),2*pi);

%% X_stim2
for i = 1:nstate
    %X_stim2{i}(:,1) = mean(coh_est)*mean((d_length))*[sin(a_unique{i}(1));cos(a_unique{i}(1))];
    X_stim2_cell{i}(:,1) = mean(coh_est,'omitnan')*0.0667*[sin(a_unique{i}(1));cos(a_unique{i}(1))];
    for j = 2:length(a_unique{i})
        %X_stim2{i}(:,j) = X_stim2{i}(:,j-1)+(1-mean(coh_est))/(n_life*ndots)*mean((d_length))*randn(2,1);
        X_stim2_cell{i}(:,j) = X_stim2_cell{i}(:,j-1)+sqrt(0.0667^2*(1-mean(coh_est,'omitnan'))/(n_life*ndots))*randn(2,1);
    end
end

X_stim2 = cell2mat(X_stim2_cell);
dir_stim2 = mod(atan2(X_stim2(1,:),X_stim2(2,:)),2*pi);


a_evidence = mod(atan2(sum_x,sum_y),2*pi);

figure
subplot(2,1,1)
plot(sum_x,'b-','LineWidth',3)
hold on
plot(X_stim2(1,:),'r-','LineWidth',3)
hold on
plot(mu(1,:),'k-','LineWidth',3)
xlabel('Frame')
ylabel('X')
set(gca,'FontSize',30)

subplot(2,1,2)
plot(sum_y,'b-','LineWidth',3)
hold on
plot(X_stim2(2,:),'r-','LineWidth',3)
hold on
plot(mu(2,:),'k-','LineWidth',3)
xlabel('Frame')
ylabel('Y')
set(gca,'FontSize',30)


figure
plot(rad2deg(unwrap(a_unique_mat)),'k-','LineWidth',2)
hold on
plot(rad2deg(unwrap(a_evidence)),'b-','LineWidth',3)
hold on
plot(rad2deg(unwrap(dir_stim2)),'r-','LineWidth',3)
legend('Veridical direction','Data stimulus','Generated stimulus')
xlabel('Frame')
ylabel('Direction')
title(['Coh = ' num2str(coh_val)])
set(gca,'FontSize',30)

%% Difference between two frames

for i = 1:length(sum_x)-1
    x_data_diff(i) = sum_x(i+1) - sum_x(i);
    y_data_diff(i) = sum_y(i+1) - sum_y(i);
    if abs(x_data_diff(i)) > 1e-3
        x_data_diff(i) = NaN;
    end
    if abs(y_data_diff(i)) > 1e-3
        y_data_diff(i) = NaN;
    end
end

for i = 1:length(X_stim2)-1
    x_stim2_diff(i) = X_stim2(1,i+1) - X_stim2(1,i);
    y_stim2_diff(i) = X_stim2(2,i+1) - X_stim2(2,i);
    if abs(x_stim2_diff(i)) > 1e-3
        x_stim2_diff(i) = NaN;
    end
    if abs(y_stim2_diff(i)) > 1e-3
        y_stim2_diff(i) = NaN;
    end
end

for i = 1:length(X_stim1)-1
    x_stim1_diff(i) = X_stim1(1,i+1) - X_stim1(1,i);
    y_stim1_diff(i) = X_stim1(2,i+1) - X_stim1(2,i);
    if abs(x_stim1_diff(i)) > 1e-3
        x_stim1_diff(i) = NaN;
    end
    if abs(y_stim1_diff(i)) > 1e-3
        y_stim1_diff(i) = NaN;
    end
end

figure
subplot(2,1,1)
plot(x_data_diff,'b-','LineWidth',3)
hold on
plot(x_stim1_diff,'r-','LineWidth',3)
legend('Data','Generated data')
xlabel('Frame')
ylabel('X diff')
set(gca,'FontSize',30)

subplot(2,1,2)
plot(y_data_diff,'b-','LineWidth',3)
hold on
plot(y_stim1_diff,'r-','LineWidth',3)
xlabel('Frame')
ylabel('Y diff')
set(gca,'FontSize',30)