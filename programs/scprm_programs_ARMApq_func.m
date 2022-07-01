% ARMA function
% Minimize the circular mean square error

function [const,js_ARMA,AIC,BIC] = scprm_programs_ARMApq_func(p,q,js_opt,a_evidence_un)

js_opt_un = js_opt;

for i = 1:p
    js_opt_un_past{i} = circshift(js_opt_un,i);
    js_opt_un_past{i}(1:max(p,q)) = [];
end

js_opt_un(1:max(p,q)) = [];


a_evidence_un_past{1} = a_evidence_un; 
a_evidence_un_past{1}(1:max(p,q)) = [];

for i = 2:q
    a_evidence_un_past{i} = circshift(a_evidence_un,i-1);
    a_evidence_un_past{i}(1:max(p,q)) = [];
end

N = length(js_opt_un_past{1});

B = @(z) 0;
for i = 1:p
    B = @(z) B(z) + z(i)*js_opt_un_past{i};
end

C = @(z) 0;
for i = 1:q
    C = @(z) C(z) + z(p+i)*a_evidence_un_past{i};
end

% Now we write it as sum of B (output terms) and C (input terms)

CMSE = @(z) mean(1 - cos(js_opt_un - B(z) - C(z)));

%z_opt = fminsearch(CMSE,rand(1,p+q));


for i = 1:10
    [z_opt{i},err_armapq(i)] = fminsearch(CMSE,rand(1,p+q));
end

i_min_armapq = find(err_armapq == min(err_armapq));
for i = 1:p+q
    const(i) = z_opt{i_min_armapq(1)}(i);
end

%{
for i = 1:p+q
    const(i) = z_opt(i);
end
%}
js_ARMA = zeros(1,length(js_opt_un));

js_ARMA(1) = js_opt_un(1);

for i = 2:p
    B_opt1 = 0;
    C_opt1 = 0;
    if q <= i
        for j = 1:q
            C_opt1 = C_opt1 + const(p+j)*a_evidence_un_past{j}(i);
        end
    elseif q > i
        for j = 1:i
            C_opt1 = C_opt1 + const(p+j)*a_evidence_un_past{j}(i);
        end
    end
    for j = 1:i-1
        B_opt1 = B_opt1 + const(j)*js_ARMA(i-j);
    end
    js_ARMA(i) = B_opt1+C_opt1;
end

% This works
for i = p+1:length(js_ARMA)
    B_opt = 0;
    C_opt = 0;
    for j = 1:q
        C_opt = C_opt + const(p+j)*a_evidence_un_past{j}(i);
    end
    for j = 1:p
        B_opt = B_opt + const(j)*js_ARMA(i-j);
    end
    js_ARMA(i) = B_opt+C_opt;
end

AIC = N*log(sum((1-cos(js_opt_un-js_ARMA)).^2)/N)+2*(p+q+1);
BIC = N*log(sum((1-cos(js_opt_un-js_ARMA)).^2)/N)+log(N)*(p+q+1);