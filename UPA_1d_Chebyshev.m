% generates Dolph-Chebyshev weights for UPA. 2N+1 elements are placed on 
%z-axis

% dayi
% 11/21/2017
% http://www.antenna-theory.com/arrays/weights/dolph.php
% http://www.antenna-theory.com/arrays/weights/dolph2.php

% close all
clear all
format compact

N = 2;
dz = 0.5; % in terms of lambda
theta = pi/180*(0:.1:180); % angle of arrival (rad)
% theta_d = 45*pi/180; % desired steering angle (rad)

%       N
% AF = SUM w_n*exp(-j*k*n*dz*cos(theta)) ; where k= 2*pi/lambda, w_n=w_-n
%      n=-N
%
%
%       N                                | w_n,n=0
%    = SUM w_n'*cos(2*n'*u) ; where w_n'=| 2*w_n,else & u=k*dz*cos(theta)/2
%     n'=0
%

syms u;
w=cell(1,N+1);
sym_keys{1}='w_0';
w{1}=sym(sym_keys{1});
AF=w{1};
for ii=1:N
    sym_keys{ii+1}=sprintf('w_%i',ii);
    w{ii+1}=sym(sym_keys{ii+1});
    AF=AF+2*w{ii+1}*cos(2*ii*u);
end
af=expand(AF);

% Substitution to have this in similar form to Chebyshev poly: cos(u)t_0=t
% t_0 is a factor to adjust SLA

SLA_dB = 30; % sidelobe attenuation
SLA_lin=20*log10(SLA_dB);
t_0=cosh(acosh(SLA_lin)/N-1);
syms t;

af=subs(af,cos(u),t/t_0);

% af=T_{2N}
T=generate_Chebyshev_sym_polynomial(2*N);
af_t=children(collect(af,t)); % this should give us N+1 coefficients
af_t=subs(af_t,t,1);
T_t=children(collect(T,t)); % equation whose terms ordered in terms of pow
T_t=subs(T_t,t,1);

for ii=1:N+1
    eqns(ii)= af_t(ii) == T_t(ii);
end

S = solve(eqns);

% now we have weights, let's substitute them in AF to see their impact
AF_num=zeros(1,length(theta));
% method-1
% for ii=0:1:N
%     AF_num=AF_num+(1+sign(ii))*double(getfield(S,sym_keys{ii+1}))...
%         *cos(2*ii*2*pi*dz*cos(theta)/2);
% end
% method-2
for ii=-1*N:1:N
    AF_num=AF_num+double(getfield(S,sym_keys{abs(ii)+1}))...
        *exp(-j*2*pi*ii*dz*cos(theta));
end

% note: both methods give the same result with expanded main lobe, which I
% am not() sure is (in)correct

% normalization of numeric AF
AF_num=AF_num/max(AF_num);
figure;plot(theta*180/pi, 20*log10(abs(AF_num(end:-1:1))))
title('Array Factor when no steering (\theta=90)')
ylabel('Gain (dB)')
xlabel('\theta')