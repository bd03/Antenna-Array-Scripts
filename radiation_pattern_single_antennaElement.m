% plots radiation pattern of a single antenna element
% according to 3GPP TS38.901 Table 7.3-1
% dayi
% 11/20/2017

clear all
close all

theta_pp=0:.1:180;
phi_pp=-180:.1:180;
A_max=30;
theta_3dB=65;
phi_3dB=65;

A_zen=-1*min(12*(theta_pp-90).^2/theta_3dB,A_max);
A_az=-1*min(12*(phi_pp/phi_3dB).^2,A_max);
A_zen1=[ones(size(theta_pp)); A_zen];
A_az1=[A_az; ones(size(phi_pp))];
A=-1*min(-1*(A_az1'*A_zen1),A_max);

At=A.';
figure;surf(phi_pp,theta_pp,At(1:length(theta_pp),1:length(phi_pp)))
shading flat
axis equal
colorbar('southoutside')
view(0,90)