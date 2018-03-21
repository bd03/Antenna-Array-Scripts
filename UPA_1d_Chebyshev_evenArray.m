% generates Dolph-Chebyshev weights for UPA. 2N elements are placed on 
%z-axis

% dayi
% 11/23/2017
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