% Compare leaky integrator with solo data.

% fname = '20211119_cla_cpr.mwk2';
% d = MW_readFile(fname, 'include', {'#stimDisplay'}, 'dotPositions');
clear all; close all;
%% At DPZ
%% On Mac, Finder, Command K, cifs://dpz.lokal/dpz
% addpath //Volumes/dpz/userinterchange/KNguyen/ 
%% On windows
% addpath Z:\userinterchange\KNguyen
% In case the network fails:
%load('//dpz.lokal/dpz/userinterchange/Felix Schneider/forKate/dots_only.mat')
%load('\\172.17.14.171\userinterchange\Felix Schneider\forKate\dots_only.mat')
addpath /Users/katenguyen/Documents/DPZ/MATLAB/Data/Solo 
addpath /Users/katenguyen/Dropbox/DPZ/MATLAB/Data
addpath /Users/katenguyen/Dropbox/DPZ/MATLAB/Programs
addpath /Users/knguyen/Dropbox/DPZ/MATLAB/Data
addpath /Users/knguyen/Dropbox/DPZ/MATLAB/Programs
addpath /Users/knguyen/Dropbox/DPZ/MATLAB/Programs/MWorks_Matlab
addpath /Users/knguyen/Dropbox/DPZ/MATLAB/Programs/MATLAB
addpath /Users/knguyen/Documents/MATLAB


load('20220214_sol_CPRsolo_block1_fxs.mat')

dt = 1/120;
n_life = 25;

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

% Convert time stamp 

for i = 1:length(dp_all_ts)-1
    js_to_fr{i} = find(and(d_jsDir_ts >= dp_all_ts(i),d_jsDir_ts < dp_all_ts(i+1)) == 1);
    js_to_fr_length(i) = length(js_to_fr{i});
end
dp_jsDir_pf{1} = cell2mat(dp_dir(js_to_fr{1}(1) - 1));
dp_jsDir_pf{2} = cell2mat(dp_dir(js_to_fr{1}(end)));
dp_jsStr_pf{1} = cell2mat(dp_str(js_to_fr{1}(1) - 1));
dp_jsStr_pf{2} = cell2mat(dp_str(js_to_fr{1}(end)));
for i = 2:length(js_to_fr)
    if length(js_to_fr{i}) == 1
        dp_jsDir_pf{i+1} = cell2mat(dp_dir(js_to_fr{i}));
        dp_jsStr_pf{i+1} = cell2mat(dp_str(js_to_fr{i}));
    %elseif length(js_to_fr{i}) == 0
    elseif isempty(js_to_fr{i}) == 1
        dp_jsDir_pf{i+1} = dp_jsDir_pf{i};
        dp_jsStr_pf{i+1} = dp_jsStr_pf{i};
    else
        %js_to_fr{i}(end)
        dp_jsDir_pf{i+1} = cell2mat(dp_dir(js_to_fr{i}(end)));
        dp_jsStr_pf{i+1} = cell2mat(dp_str(js_to_fr{i}(end)));
        %dp_jsDir_pf{i+1} = dp_jsDir_pf{i}(end);
        %dp_jsStr_pf{i+1} = dp_jsStr_pf{i}(end);
    end
end

%bl = 6;
bl = 46; % block no.
coh_val = dc{bl};
idx.bl = and(d_rdp_dir_ts>=dc_ts(bl),d_rdp_dir_ts< dc_ts(bl+1)); % indices of steady states within that block
ss = find(idx.bl == 1); % indices steady states in that block
nstate = length(ss);
substate = and(dp_all_ts>=d_rdp_dir_ts(ss(1)),dp_all_ts<d_rdp_dir_ts(ss(end)+1));
ss_ts = find(substate == 1);
dp_ts = dp_all_ts(ss_ts);

dp = dp_all(ss_ts);
js_dir_ss = dp_jsDir_pf(ss_ts);
js_str_ss = dp_jsStr_pf(ss_ts);

%% Remove the last point, need to think about it
js_dir_ss(1) = [];
js_str_ss(1) = [];

js_dir_mat = deg2rad(cell2mat(js_dir_ss));

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

dt = mean(mode(d_length))/10;

a_true = zeros(length(dir),length(dir{1}));
for i = 1:length(dir) % direction
    for j = 1:length(dir{1}) % frame
        a_true(i,j) = dir{i}(j);
    end
end

dir_test = dir;

for i = 1:nstate
    st = and(dp_all_ts>=d_rdp_dir_ts(ss(i)),dp_all_ts<d_rdp_dir_ts(ss(i)+1));
    t_index{i} = find(st == 1); 
    a_unique{i} = deg2rad(cell2mat(rdp_dir(ss(i))))*ones(1,length(t_index{i}));
end

a_unique_mat = cell2mat(a_unique);

% Remove dots that disappear from d_length
for i = 1:ndots
    for j = 1:iFrame-1
        if abs(d_length(i,j)-mode(d_length(:,j))) > 4.9e-4
            dir{i}(j) = NaN;
            x_resultant{i}(j) = 0;
            y_resultant{i}(j) = 0;
            a_true(i,j) = NaN;
        end
    end
end

for i = 1:ndots
    x_speed(i,:) = x_resultant{i};
    y_speed(i,:) = y_resultant{i};
end

for j = 1:iFrame-1
    ndots_avail(j) = ndots - length(find(x_speed(:,j) == 0));
end

sum_x = sum(x_speed)./ndots_avail;% x-evidence per frame;
sum_y = sum(y_speed)./ndots_avail;% y-evidence per frame;

for i = 1:length(sum_x)
    dir_evidence(i) = mod(atan2(sum_x(i),sum_y(i)),2*pi);
end

a_evidence = mod(atan2(sum_x,sum_y),2*pi);

%% Old Stimulus

for i = 1:length(a_unique_mat)-1
    mu(:,i) = coh_val*mean(mode(d_length))*[sin(a_unique_mat(i));cos(a_unique_mat(i))];
end
VX = 0.5*(mean(mode(d_length)))^2*(1-coh_val)/ndots;
S = [VX 0;0 VX];

for i = 1:length(mu)
    X_stim(:,i) = mvnrnd(mu(:,i),S);
end

dir_stim = mod(atan2(X_stim(1,:),X_stim(2,:)),2*pi);


%% X_stim2 that includes the past stimulus

%% X_stim2
for i = 1:nstate
    %X_stim2{i}(:,1) = mean(coh_est)*mean((d_length))*[sin(a_unique{i}(1));cos(a_unique{i}(1))];
    X_stim2_cell{i}(:,1) = coh_val*0.0667*[sin(a_unique{i}(1));cos(a_unique{i}(1))];
    for j = 2:length(a_unique{i})
        %X_stim2{i}(:,j) = X_stim2{i}(:,j-1)+(1-mean(coh_est))/(n_life*ndots)*mean((d_length))*randn(2,1);
        X_stim2_cell{i}(:,j) = X_stim2_cell{i}(:,j-1)+sqrt(0.0667^2*(1-coh_val)/(n_life*ndots))*randn(2,1);
    end
end

X_stim2 = cell2mat(X_stim2_cell);
dir_stim2 = mod(atan2(X_stim2(1,:),X_stim2(2,:)),2*pi);
dir_stim2(end) = [];

figure(1)
plot(rad2deg(unwrap(dir_stim)),'g-','LineWidth',3)
hold on
plot(rad2deg(unwrap(dir_stim2)),'r-','LineWidth',3)
hold on
plot(rad2deg(unwrap(a_evidence)),'b-','LineWidth',3)
hold on
plot(rad2deg(unwrap(a_unique_mat)),'k-','LineWidth',3)
xlabel('Frame')
ylabel('Direction')
legend('Stim 1','Stim 2','Data','Veridical')
title(['Coh = ' num2str(coh_val)])
set(gca,'FontSize',30)


% Minimize distance between direction evidence and joystick direction
%f = @(TD,x) mean(1 - cos(a_unique_mat(1:end-1)-circshift(js_dir_mat,-round(TD/dt))));
f = @(TD,x) mean(1 - cos(x-circshift(js_dir_mat,-round(TD/dt))));
f_dir = @(TD) f(TD,dir_evidence);
f_stim = @(TD) f(TD,dir_stim2);

for i = 1:10
    %[TD_min1(i),min_dir(i)] = fminsearchbnd(f_dir,rand,0);
    [TD_min2(i),min_stim(i)] = fminsearchbnd(f_stim,rand,0);
end

i_min_stim = find(min_stim == min(min_stim));

TD_stim = TD_min2(i_min_stim(1));


%js_dir_opt1 = circshift(js_dir_mat,-round(TD_dir/dt));
js_dir_opt2 = circshift(js_dir_mat,-round(TD_stim/dt));


x_js_stim = sin(((js_dir_opt2)));
y_js_stim = cos(((js_dir_opt2)));

a_unique_un = unwrap(a_unique_mat(1:end-round(TD_stim/dt)));
js_dir_un = unwrap(js_dir_opt2(1:end-round(TD_stim/dt)));

for i = 1:length(js_dir_un)
    if js_dir_un(i) - a_unique_un(i) > pi
        js_dir_un(i) = js_dir_un(i) - 2*pi;
    elseif js_dir_un(i) - a_unique_un(i) < -pi
        js_dir_un(i) = js_dir_un(i) + 2*pi;
    end
end


figure
plot(rad2deg(a_unique_un),'k-','LineWidth',3)
hold on
plot(rad2deg(unwrap(dir_evidence(1:end-round(TD_stim/dt)))),'g-','LineWidth',3)
hold on
plot(rad2deg(unwrap(dir_stim2(1:end-round(TD_stim/dt)))),'r-','LineWidth',3)
hold on
plot(rad2deg(js_dir_un),'b-','LineWidth',3)
xlabel('Frame')
ylabel('Direction')
legend('Veridical dir.','Stim. from data','Generated stim.','Joystick direction')
title(['Coh = ' num2str(coh_val)])
set(gca,'FontSize',30)

x_js_stim(end-round(TD_stim/dt)+1:end) = [];
y_js_stim(end-round(TD_stim/dt)+1:end) = [];
X_stim2(:,end-round(TD_stim/dt)+1:end) = [];
sum_x(end-round(TD_stim/dt)+1:end) = [];
sum_y(end-round(TD_stim/dt)+1:end) = [];
X_stim2(:,end) = [];


%% Leaky integrator using x_sum vs. X_stim
% Optimaize between joystick and X_stim first, then use that optimal values
% to find goodness-of-fit between joystick and x_sum
% function X_leak = scprm_programs_leak_func(X_stim,ti,dt,lambda,T_delay)

%MSE = @(lambda,TD,x,y,xj,yj,k) mean((norm([xj;yj] - k*leak_position_func(x,y,lambda,TD,dt))).^2);

ti = 0:dt:(length(x_js_stim)-1)*dt;
MSE = @(lambda,TD,x,y,xj,yj,k) mean((norm([xj;yj] - k*scprm_programs_leak_func([x;y],ti,dt,lambda,TD))).^2);


MSE_Xstim = @(z) MSE(z(1),z(2),X_stim2(1,:),X_stim2(2,:),x_js_stim,y_js_stim,z(3));

z_opt2 = fminsearchbnd(MSE_Xstim,rand(1,3),[-100,0,-100]);
l_opt_Xstim = z_opt2(1);
TD_opt_Xstim = z_opt2(2);
k_opt_Xstim = z_opt2(3);


X_leak_Xstim_opt = k_opt_Xstim*scprm_programs_leak_func(X_stim2,ti,dt,l_opt_Xstim,TD_opt_Xstim);
X_leak_sumx_opt = k_opt_Xstim*scprm_programs_leak_func([sum_x;sum_y],ti,dt,l_opt_Xstim,TD_opt_Xstim);



figure(2)
subplot(2,1,1)
plot(x_js_stim,'k-','LineWidth',3)
hold on
plot(X_leak_Xstim_opt(1,:),'b-','LineWidth',3)
xlabel('Frame')
ylabel('X')
set(gca,'FontSize',30)

subplot(2,1,2)
plot(y_js_stim,'k-','LineWidth',3)
hold on
plot(X_leak_Xstim_opt(2,:),'b-','LineWidth',3)
xlabel('Frame')
ylabel('Y')
set(gca,'FontSize',30)

figure(3)
subplot(2,1,1)
plot(x_js_stim,'k-','LineWidth',3)
hold on
plot(X_leak_sumx_opt(1,:),'b-','LineWidth',3)
xlabel('Frame')
ylabel('X')
set(gca,'FontSize',30)

subplot(2,1,2)
plot(y_js_stim,'k-','LineWidth',3)
hold on
plot(X_leak_sumx_opt(2,:),'b-','LineWidth',3)
xlabel('Frame')
ylabel('Y')
set(gca,'FontSize',30)

%------------------------------------

a_leak_dir = mod(atan2(X_leak_sumx_opt(1,:),X_leak_sumx_opt(2,:)),2*pi);
a_leak_stim = mod(atan2(X_leak_Xstim_opt(1,:),X_leak_Xstim_opt(2,:)),2*pi);

a_leak_dir_un = unwrap(a_leak_dir);
a_leak_stim_un = unwrap(a_leak_stim);

js_dir_opt2(end-round(TD_stim/dt)+1:end) = [];

js_dir_un2 = unwrap(js_dir_opt2);

a_unique_un = unwrap(a_unique_mat);
a_unique_un(end-round(TD_stim/dt)+1:end) = [];

for i = 1:length(js_dir_un2)
    if js_dir_un2(i) - a_unique_un(i) > pi
        js_dir_un2(i) = js_dir_un2(i) - 2*pi;
    elseif js_dir_un2(i) - a_unique_un(i) < -pi
        js_dir_un2(i) = js_dir_un2(i) + 2*pi;
    end
end

for i = 1:length(js_dir_un2)
    if a_leak_dir_un(i) - a_unique_un(i) > pi
        a_leak_dir_un(i) = a_leak_dir_un(i) - 2*pi;
    elseif a_leak_dir_un(i) - a_unique_un(i) < -pi
        a_leak_dir_un(i) = a_leak_dir_un(i) + 2*pi;
    end
end

for i = 1:length(a_leak_stim_un)
    if a_leak_stim_un(i) - a_unique_un(i) > pi
        a_leak_stim_un(i) = a_leak_stim_un(i) - 2*pi;
    elseif a_leak_stim_un(i) - a_unique_un(i) < -pi
        a_leak_stim_un(i) = a_leak_stim_un(i) + 2*pi;
    end
end

%a_leak_stim(2191:2193) = [];
a_leak_stim_un(1) = [];
a_leak_dir_un(1) = [];
figure(4)
%subplot(2,1,1)
plot(rad2deg(js_dir_un2),'g-','LineWidth',3)
hold on
plot(rad2deg(a_leak_dir_un),'b-','LineWidth',3)
hold on
plot(rad2deg(a_leak_stim_un),'r-','LineWidth',3)
hold on
plot(rad2deg(a_unique_un),'k-','LineWidth',3)
xlabel('Frame')
ylabel('Direction')
legend('JS direction','Frame-by-frame data','Generated data','Veridical')
title(['Coh = ' num2str(coh_val)])
set(gca,'FontSize',30)

CMSE_dir = mean(1-cos(js_dir_opt2-a_leak_dir));
CMSE_stim = mean(1-cos(js_dir_opt2-a_leak_stim));

%% Compute AIC 

N = length(js_dir_opt2);
AIC_dir = N*log(CMSE_dir/N) + 2*(2+1);
AIC_stim = N*log(CMSE_stim/N) + 2*(2+1);


figure(5)
%subplot(2,1,1)
plot(rad2deg(js_dir_un2),'k-','LineWidth',3)
hold on
plot(rad2deg(a_leak_dir_un),'b-','LineWidth',3)
hold on
plot(rad2deg(a_leak_stim_un),'r-','LineWidth',3)
xlabel('Frame')
ylabel('Direction')
legend('JS direction','Frame-by-frame data','Generated data')
title(['Coh = ' num2str(coh_val)])
set(gca,'FontSize',30)
