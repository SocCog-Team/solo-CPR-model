% X_sig = (x,y)     : veridical direction
% X_stim = (x',y')  : stimulus direction (veridical + noise)
% ti                : time
% dt = 1/Frame rate 
% T_delay           : perceptual + motor delay
% lambda            : discount rate

function X_leak = scprm_programs_leak_func(X_stim,ti,dt,lambda,T_delay)

% Each frame lasts for 10ms, perceptual delay is 100ms ~ 10 frames
i_delay = round(T_delay/dt);

X_leak = zeros(2,length(ti));

for i = 2:length(ti)-(i_delay)
    X_leak(:,i+i_delay) = X_leak(:,i+i_delay-1)+X_stim(:,i)*dt-lambda*X_leak(:,i+i_delay-1)*dt;
end
    
X_leak = lambda.*X_leak;