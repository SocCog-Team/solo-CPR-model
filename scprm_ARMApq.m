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
addpath /Users/katenguyen/Documents/DPZ/MATLAB/Data/Solo % DPZ windows
addpath /Users/katenguyen/Dropbox/DPZ/MATLAB/Data
addpath /Users/katenguyen/Dropbox/DPZ/MATLAB/Programs
addpath /Users/knguyen/Dropbox/DPZ/MATLAB/Data
addpath /Users/knguyen/Dropbox/DPZ/MATLAB/Programs
addpath /Users/knguyen/Dropbox/DPZ/MATLAB/Programs/MWorks_Matlab
addpath /Users/knguyen/Dropbox/DPZ/MATLAB/Programs/MATLAB
addpath /Users/knguyen/Documents/MATLAB


load('20220214_sol_CPRsolo_block1_fxs.mat')

dt = 1/120;

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
    elseif isempty(js_to_fr{i}) == 1
        dp_jsDir_pf{i+1} = dp_jsDir_pf{i};
        dp_jsStr_pf{i+1} = dp_jsStr_pf{i};
    else
        dp_jsDir_pf{i+1} = cell2mat(dp_dir(js_to_fr{i}(end)));
        dp_jsStr_pf{i+1} = cell2mat(dp_str(js_to_fr{i}(end)));
    end
end

%bl = 6;
bl = 2; % block no.
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

js_dir_mat = unwrap(deg2rad(cell2mat(js_dir_ss)));

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

for i = 1:length(dir_evidence)
    if (dir_evidence(i) - js_dir_mat(i)) > pi + 0.1
        js_dir_mat(i) = js_dir_mat(i) + 2*pi;
    elseif dir_evidence(i) - js_dir_mat(i) < -pi - 0.1
        js_dir_mat(i) = js_dir_mat(i) - 2*pi;
    end
end

a_unique_mat(end) = [];

% Minimize distance between direction evidence and joystick direction
f = @(TD) mean(1 - cos(a_unique_mat-circshift(js_dir_mat,-round(TD/dt))));
TD_opt = fminsearchbnd(f,rand,0);

js_opt = circshift(js_dir_mat,-round(TD_opt/dt));

js_opt(end-round(TD_opt/dt)+1:end) = [];
a_unique_mat(end-round(TD_opt/dt)+1:end) = [];
dir_evidence(end-round(TD_opt/dt)+1:end) = [];
sum_x(end-round(TD_opt/dt)+1:end) = [];
sum_y(end-round(TD_opt/dt)+1:end) = [];

a_unique_un = unwrap(a_unique_mat);
a_evidence_un = unwrap(dir_evidence);
js_opt_un = unwrap(js_opt);

for i = 1:length(js_opt_un)
    if js_opt_un(i) - a_unique_un(i) > pi
        js_opt_un(i) = js_opt_un(i) - 2*pi;
    elseif js_opt_un(i) - a_unique_un(i) < -pi
        js_opt_un(i) = js_opt_un(i) + 2*pi;
    end
end

figure(1)
plot(rad2deg(a_evidence_un),'r-','LineWidth',3)
hold on
plot(rad2deg(a_unique_un),'k-','LineWidth',3)
hold on
plot(rad2deg(js_opt_un),'b-','LineWidth',3)
legend('Stimulus','Veridical','Joystick direction')
title(['Coh = ' num2str(coh_val)])
xlabel('Frame')
ylabel('Direction')
set(gca,'FontSize',30)

%-----------------------------
% Leaky integrator
%-----------------------------
% scprm_programs_leak_func(X_stim,ti,dt,lambda,T_delay)
% function X_leak = leak_position_func(x,y,lambda,TD,dt)

ti = 0:dt:dt*(length(sum_x)-1);
MSE = @(lambda,TD,x,y,xj,yj,k) mean((norm(k*[xj;yj] - scprm_programs_leak_func([x;y],ti,dt,lambda,TD))).^2);

x_js_dir = sin(js_opt_un);
y_js_dir = cos(js_opt_un);

MSE_sumx = @(m) MSE(m(1),0,sum_x,sum_y,x_js_dir,y_js_dir,m(2));
m_opt = fminsearchbnd(MSE_sumx,rand(1,2),zeros(1,2));

l_opt = m_opt(1);
k_opt = m_opt(2);

X_leak = scprm_programs_leak_func([sum_x;sum_y],ti,dt,l_opt,0);

a_leak = mod(atan2(X_leak(1,:),X_leak(2,:)),2*pi);

figure(2)
subplot(2,1,1)
plot(sum_x,'k-','LineWidth',3)
hold on
plot(k_opt*x_js_dir,'r-','LineWidth',3)
hold on
plot(X_leak(1,:),'b-','LineWidth',3)
legend('Stimulus','Joystick direction','Leaky integrator')
xlabel('Frame')
ylabel('X')
title(['Coh = ' num2str(coh_val)])
set(gca,'FontSize',30)

subplot(2,1,2)
plot(sum_y,'k-','LineWidth',3)
hold on
plot(k_opt*y_js_dir,'r-','LineWidth',3)
hold on
plot(X_leak(2,:),'b-','LineWidth',3)
xlabel('Frame')
ylabel('X')
set(gca,'FontSize',30)

a_leak_un = unwrap(a_leak);

for i = 1:length(a_leak_un)
    if a_leak_un(i) - a_unique_un(i) > pi
        a_leak_un(i) = a_leak_un(i) - 2*pi;
    elseif a_leak_un(i) - a_unique_un(i) < -pi
        a_leak_un(i) = a_leak_un(i) + 2*pi;
    end
end

figure(3)
plot(rad2deg(a_unique_un),'k-','LineWidth',3)
hold on
plot(rad2deg(unwrap(js_opt)),'b-','LineWidth',3)
hold on
plot(rad2deg(a_leak_un(2:end)),'r-','LineWidth',3)
legend('Veridical','Joystick','Leaky integrator')
title(['Coh = ' num2str(coh_val)])
xlabel('Frame')
ylabel('Direction')
set(gca,'FontSize',30)

%-----------------------------
% ARMA(p,q)
% output = b1*output(-1)+b2*output(-2)+ c1*input+c2*input(-1)
%-----------------------------
p = 2; q = 2;

[const,js_ARMA,AIC,BIC] = scprm_programs_ARMApq_func(p,q,js_opt_un,a_evidence_un);

js_ARMA_un = unwrap(js_ARMA);

for i = 1:length(js_ARMA_un)
    if js_ARMA_un(i) - a_unique_un(i) > 2*pi
        js_ARMA_un(i) = js_ARMA_un(i) - 2*pi;
    elseif js_ARMA_un(i) - a_unique_un(i) < -2*pi
        js_ARMA_un(i) = js_ARMA_un(i) + 2*pi;
    end
end

figure(4)
plot(rad2deg(js_opt_un(max(p,q)+1:end)),'b-','LineWidth',3)
hold on
plot(rad2deg(js_ARMA_un),'r-','LineWidth',3)
legend('Joystick direction','ARMA(2,2)')
title(['Coh = ' num2str(coh_val)])
xlabel('Frame')
ylabel('Direction')
set(gca,'FontSize',30)

%-----------------------------
% ARMA(1,1) for leaky integrator
% output = b*output(-1) + c*input
%-----------------------------

a_leak_un_1 = circshift(a_leak_un,1);
a_leak_un_1(1) = [];


[const_leak,js_ARMA_leak,AIC_leak,BIC_leak] = scprm_programs_ARMApq_func(1,1,a_leak_un,a_evidence_un);


js_ARMA_leak_un = unwrap(js_ARMA_leak);
for i = 1:length(a_leak_un)
    if js_ARMA_leak_un(i) - a_unique_un(i) > pi
        js_ARMA_leak_un(i) = js_ARMA_leak_un(i) - 2*pi;
    elseif js_ARMA_leak_un(i) - a_unique_un(i) < -pi
        js_ARMA_leak_un(i) = js_ARMA_leak_un(i) + 2*pi;
    end
end

figure(5)
plot(rad2deg(js_opt_un),'k-','LineWidth',3)
hold on
plot(rad2deg(a_leak_un(2:end)),'b-','LineWidth',3)
hold on
legend('Joystick','Leaky integrator','ARMA(1,1)')
title(['Coh = ' num2str(coh_val)])
xlabel('Frame')
ylabel('Direction')
set(gca,'FontSize',30)