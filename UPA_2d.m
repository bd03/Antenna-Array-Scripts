% calculates the output of the isotropic receiver antennas placed uniformly
% over xy-plane

% dayi
% 11/21/2017
% http://www.antenna-theory.com/arrays/weights/twoDuniform.php

close all
clear all

% 2d ULA assuming that (2M+1)X(2N+1) AEs are placed at (x_mn,y_mn,0),
% m=1,...,2M+1 and n=1,...,2N+1 where AE with index (M,N) is placed at 
% (0, 0, 0).

M = 1; % number of rows (on y-axis)
N = 1; % number of columns (x-axis)
dx = 0.5; % in terms of lambda
dy = 0.7; % in terms of lambda
theta = pi/180*(0:.1:180); % angle of arrival (rad)
theta_d = 45*pi/180; % desired steering angle (rad)
phi=pi/180*(-180:.1:180);
phi_d=45*pi/180;
pseudo_lambda=1;

% wave vector (steering)
k_d=2*pi/pseudo_lambda*[sin(theta_d)*cos(phi_d); sin(theta_d)*sin(phi_d); cos(theta_d)];
% AE locations
r=[repmat((-1*N:1:N)*dx,1,2*M+1);...
    reshape(repmat((-1*N:1:N).'*dy, 1,2*N+1).', 1,(2*M+1)*(2*N+1));...
    zeros(1,(2*M+1)*(2*N+1))];

w_temp=(exp(j*k_d.'*r));
w=reshape(w_temp,2*M+1,2*N+1).';

% Array factor computation
AF=zeros(length(theta),length(phi));
% AE locations reshaped
r=reshape(r,3,2*M+1,2*N+1);
% All wave vectors
k=2*pi/pseudo_lambda*[reshape(sin(theta).'*cos(phi),1,length(theta)*length(phi));...
    reshape(sin(theta).'*sin(phi),1,length(theta)*length(phi));...
    repmat(cos(theta),1,length(phi))];
% steering vector v(k)=[exp(-j*k*r_n)]
% AF=weight_vector*steering_vector
for mm=1:(2*M+1)
    for nn=1:(2*N+1)
        v(mm,nn,:)=exp(-j*k.'*r(:,mm,nn));
        AF=AF+w(mm,nn)*reshape(v(mm,nn,:),length(theta),length(phi));
    end
end
AF=AF/max(max(AF));
AF_dB=20*log10(abs(AF));

% just cosmetics
indices = find(AF<-30);
AF(indices)=-30;

figure;surf(phi,theta,AF_dB(1:length(theta),1:length(phi)))
shading flat
axis equal
colorbar('southoutside')
view(0,90)
title('Array Factor')
xlabel('\theta')