% Toy model using scprm_programs_leak_func.m
% X_leak = scprm_programs_leak_func(X_stim,ti,dt,lambda,T_delay)
clear all; close all;

dt = 1/120;
ti = 0:dt:2; % 2 seconds
a_veridical = 75; % in degrees
coh = 0.8;
ndots = 500;

x_veridical = sind(a_veridical);
y_veridical = cosd(a_veridical);


mu = [x_veridical; y_veridical];

X_sig = mu*ones(1,length(ti));
VX = 0.5*(1-coh)/ndots;
S = [VX 0;0 VX]; % covariance matrix

for i = 1:length(ti)
    X_stim(:,i) = mvnrnd(mu,S);
end

%X_leak = @(lambda,T_delay) scprm_programs_leak_func(X_stim,ti,dt,lambda,T_delay);

err = @(C) scprm_programs_MSEfunc(X_sig,X_stim,C(1),ti,dt,C(2),C(3));

C_opt = fminsearchbnd(err,rand(1,3),zeros(1,3));

k_opt = C_opt(1);
l_opt = C_opt(2);
TD_opt = C_opt(3);

X_leak_opt = k_opt*scprm_programs_leak_func(X_stim,ti,dt,l_opt,TD_opt);

figure
subplot(2,1,1)
plot(ti,X_sig(1,:),'k-','LineWidth',3)
hold on
plot(ti,X_leak_opt(1,:),'b-','LineWidth',3)
xlabel('Time')
ylabel('X')
set(gca,'FontSize',30)

subplot(2,1,2)
plot(ti,X_sig(2,:),'k-','LineWidth',3)
hold on
plot(ti,X_leak_opt(2,:),'b-','LineWidth',3)
xlabel('Time')
ylabel('Y')
set(gca,'FontSize',30)

a_leak = rad2deg(mod(atan2(X_leak_opt(1,:),X_leak_opt(2,:)),2*pi));

figure
plot(ti,a_veridical*ones(1,length(ti)),'k-','LineWidth',3)
hold on
plot(ti,a_leak,'b-','LineWidth',3)
xlabel('Time')
ylabel('Direction')
set(gca,'FontSize',30)