% Compute optimal leak constant lambda
% The higher the lambda, the more information the subject discounts

clear all; close all;

addpath /Users/katenguyen/Documents/DPZ/MATLAB/Data/Solo % path to solo data
addpath /Users/katenguyen/Dropbox/DPZ/MATLAB/Programs % path to program

ndots = 100;
dt = 0.010;
%coh = 0.8;
coh = 0.25;
%coh = 0.1:0.1:0.8;
nsim = 100;
T_switch = 1;
%T_switch = linspace(1,2.5,20);
%T_switch = 1:0.2:2.5;
n_switch = 100;


for kk = 1:length(T_switch)
for k = 1:length(coh)

VX = 0.5*(1-coh(k))/ndots;
S = [VX 0; 0 VX];

T = n_switch*T_switch(kk);
ts = 0:dt:T_switch(kk);
ti = 0:dt:T;


for i = 1:nsim
    X_stim = zeros(2,length(ti));
    X_sig = []; 
    for j = 1:n_switch
        a = 2*pi*rand;
        X_sig = [X_sig repmat([coh(k)*cos(a);coh(k)*sin(a)],1,length(ts)-1)];
    end
    X_sig = [X_sig [coh(k)*cos(a);coh(k)*sin(a)]];
    for i_stim = 1:length(ti)
        X_stim(:,i_stim) = mvnrnd(X_sig(:,i_stim),S)';
    end
    [l(i),MSE(i)] = fminbnd(@(lambda)MSE_func(X_sig,X_stim,ti,dt,lambda,0),0.1,100);
end

l_mean(k,kk) = mean(l);
l_std(k,kk) = std(l);
MSE_mean(k,kk) = mean(MSE);

end
end

T = 0:dt:5;
t_switch = 0:dt:T_switch;
a = 2*pi*rand(1,5);
for i = 1:5
    mu(:,i) = [coh*cos(a(i));coh*sin(a(i))];
end
X_stim1(:,1) = mvnrnd(mu(:,1),S)';
X_sig1(:,1) = mu(:,1);
for j = 1:5
    for i = (j-1)*(length(t_switch)-1)+2:j*(length(t_switch)-1)+1
        X_stim1(:,i) = mvnrnd(mu(:,j),S)';
        X_sig1(:,i) = mu(:,j);
    end
end

X_leak1 = scprm_programs_leak_func(X_stim1,T,dt,l_mean,0);

figure(1)
subplot(2,1,1)
title(['Coh = ' num2str(coh)])
plot(T,X_stim1(1,:),'LineWidth',3,'Color','r')
hold on
plot(T,X_sig1(1,:),'LineWidth',3,'Color','k')
hold on
plot(T,X_leak1(1,:),'LineWidth',3,'Color','b')
xlabel('Time')
ylabel('X')
set(gca,'FontSize',25)

subplot(2,1,2)
title(['Coh = ' num2str(coh)])
plot(T,X_stim1(2,:),'LineWidth',3,'Color','r')
hold on
plot(T,X_sig1(2,:),'LineWidth',3,'Color','k')
hold on
plot(T,X_leak1(2,:),'LineWidth',3,'Color','b')
xlabel('Time')
ylabel('Y')
set(gca,'FontSize',25)

figure(2)
s1 = subplot(1,2,1);
set(gca,'FontSize',20)
s2 = subplot(1,2,2);
set(gca,'FontSize',20)
xlim(s1,[-2*coh 2*coh])
ylim(s1,[-2*coh 2*coh])
xlim(s2,[-2*coh 2*coh])
ylim(s2,[-2*coh 2*coh])
xlabel(s1,'X')
xlabel(s2,'X')
ylabel(s1,'Y')
ylabel(s2,'Y')
hold(s1,'on')
hold(s2,'on')

figure(2)
for i = 1:length(T)
    plot(s1,X_sig1(1,i),X_sig1(2,i),'k.','MarkerSize',30)
    plot(s1,X_stim1(1,i),X_stim1(2,i),'r.','MarkerSize',10)
    plot(s2,X_leak1(1,i),X_leak1(2,i),'b.','MarkerSize',10)
    plot(s2,X_sig1(1,i),X_sig1(2,i),'k.','MarkerSize',30)
    pause(dt)
    %F(i) = getframe(gcf);
end
%video = VideoWriter('coh025','MPEG-4');
%open(video)
%writeVideo(video,F)
%close(video)

figure(3)
subplot(1,2,1)
plot(X_stim1(1,:),X_stim1(2,:),'Color','r','LineWidth',1)
hold on
plot(X_stim1(1,:),X_stim1(2,:),'r.','MarkerSize',10)
hold on
plot(X_sig1(1,:),X_sig1(2,:),'k.','MarkerSize',30)
xlim([-2*coh 2*coh])
ylim([-2*coh 2*coh])
xlabel('X')
ylabel('Y')
title(['Coh = ', num2str(coh)])
set(gca,'FontSize',20)

subplot(1,2,2)
plot(X_leak1(1,:),X_leak1(2,:),'Color','b','LineWidth',1)
hold on
plot(X_leak1(1,:),X_leak1(2,:),'b.','MarkerSize',10)
hold on
plot(X_sig1(1,:),X_sig1(2,:),'k.','MarkerSize',30)
xlim([-2*coh 2*coh])
ylim([-2*coh 2*coh])
xlabel('X')
ylabel('Y')
title(['Coh = ', num2str(coh)])
set(gca,'FontSize',20)

function err = MSE_func(X_sig,X_stim,ti,dt,lambda,T_delay)

X_leak = scprm_programs_leak_func(X_stim,ti,dt,lambda,T_delay);
err = mean((norm(X_sig-X_leak)).^2);

end