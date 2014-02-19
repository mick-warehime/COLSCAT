close all; clear all; clc

% scattering dynamics of H+H2
mh = 1836.15264;          % mass of hydrogen in au
m  = [mh,mh,mh];          % [H,H,H]
E = 0.01:0.0001:0.06;     % energy range
initState = 0;            % initial state
 
% calculate mesh with approximately 2000 points
% mesh = hh2_mesh(2000,1);

% or load premade mesh with about 2000 points
load hh2_4218

% determine scattering probabilities and wave functions
[R,T,psi,mesh,jac,viba,vibc] = colscat(m,E,initState,1,mesh,1,'vhh2');

% collision energy in eV
ecol = (E-viba{1}.e(1))*27.211;

% plot reflection coefficients for first three vibrational states
figure(1)
plot(ecol,R(:,1:3))
legend('0-0','0-1','0-2')
title('Reflection Probabilities for H + H_2(v=0)','fontsize',24)
set(gca,'fontsize',20)
axis([ecol([1 end]) 0 1])
set(gca,'xtick',0:.25:1.5)
xlabel('Collision Energy, eV')

% plot reflection coefficients for first three vibrational states
figure(2)
plot(ecol,T(:,1:3))
legend('0-0','0-1','0-2')
title('Transmission Probabilities for H + H_2(v=0)','fontsize',24)
set(gca,'fontsize',20)
axis([ecol([1 end]) 0 1])
set(gca,'xtick',0:.25:1.5)
xlabel('Collision Energy, eV')
