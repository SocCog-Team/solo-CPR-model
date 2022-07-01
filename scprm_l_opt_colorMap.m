% Optimize lambda using fminbnd
% Use leak_func to get output of leaky integrator
% [X_stim,X_leak] = leak_func(coh,ndots,T,T_switch,dt,lambda,alpha)
% MSE_func(X_sig,coh,ndots,T,T_switch,dt,lambda,alpha)

clear all; close all;

ndots = 100;
dt = 0.010;
coh = 0.1:0.1:0.9;
nsim = 100;
T_switch = 1:0.2:2.7;
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
    l(i) = fminbnd(@(lambda)MSE_func(X_sig,X_stim,ti,dt,lambda,0),0.1,100);
end

l_mean(k,kk) = mean(l);
l_std(k,kk) = std(l);

end
end

figure(1), pcolor(coh,T_switch,l_mean'), shading flat, colormap(hot);
colorbar
xlabel('Coherence')%,'fontsize',30,'interpreter','latex')
ylabel('Switch rate')%,'fontsize',30,'interpreter','latex')
set(gca,'FontSize',30)

function err = MSE_func(X_sig,X_stim,ti,dt,lambda,T_delay)

X_leak = scprm_programs_leak_func(X_stim,ti,dt,lambda,T_delay);
err = mean((norm(X_sig-X_leak)).^2);

end