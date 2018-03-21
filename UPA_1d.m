% calculates the output of the isotropic receiver antennas placed uniformly
% on z-axis

% dayi
% 11/20/2017
% http://www.antenna-theory.com/arrays/main.php

close all
clear all

%% output of isotropic receiver antennas
% 1d ULA assuming that 2N+1 AEs are placed at (0,0,z_n), n=1,...,2N+1 where
% AE with index N is placed at (0, 0, 0).

N = 2;
dz = 0.5; % in terms of lambda
theta = pi/180*(0:.1:180); % angle of arrival (rad)
theta_d = 45*pi/180; % desired steering angle (rad)

% contribution from each antenna ii: (i=1,...,2N-1)
% Xi=exp{-j*2*pi/lambda*cos(theta)*(-N*dz + ii*dz))}

Y=zeros(size(theta));
w_ii=zeros(0,0);
for ii=-1*N:1:N
    Xii=exp(-j*2*pi*cos(theta)*(ii*dz));
    Y = Y + Xii;
%     w_ii(end+1,:)=Xii;
    w_ii(end+1)=exp(-j*2*pi*cos(theta_d)*(ii*dz)); % weights for desired steering
end
% normalization of Y
Y=Y/max(Y);

figure;plot(theta*180/pi, 20*log10(abs(Y)))
title('Array Factor when no steering (\theta=90)')
ylabel('Gain (dB)')
xlabel('\theta')

%% Array factor computation
% http://www.antenna-theory.com/arrays/weights/main.php

phi=0;
pseudo_lambda=1;
% wave vector
k=2*pi/pseudo_lambda*[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
% AE locations
r=[zeros(1,2*N+1);zeros(1,2*N+1);(-1*N:1:N)*dz*pseudo_lambda];
% steering vector v(k)=[exp(-j*k*r_n)]
% AF=weight_vector*steering_vector
AF=zeros(1,length(theta));
for ii=1:(2*N+1)
    v(ii,:)=exp(-j*k.'*r(:,ii));
    AF=AF+w_ii(ii)*v(ii,:);
end
% normalization of AF
AF=AF/max(AF);
figure;plot(theta*180/pi, 20*log10(abs(AF(end:-1:1))))
title('Array Factor with steering to \theta=\pi/4')
ylabel('Gain (dB)')
xlabel('\theta')