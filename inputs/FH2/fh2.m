close all; clear all; clc

% reactive scattering of F+H2
mh = 1836.15264;             % mass of hydrogen in au
m  = mh*[19/1.008 1 1];      % [F,H,H]
E = 0.0099:0.0001:0.0282;    % energy range
state = 0;                   % initial state
vinput = 'fh2_muck';         % potential energy function

% calculate new mesh
% mesh = fh2_mesh(2000,1);

% or load existing mesh
load fh2_8301

%% Call main Scattering Routine
[R,T,psi,mesh,jac,viba,vibc,v] = colscat(m,E,state,1,mesh,1,vinput);
 
 % collision energy in eV
ecol = (E-viba{1}.e(1))*27.211;
 
ll = 73; % pick a scattering energy to analyze in figures 2-4
 
figure(1);
hold all;
plot(ecol,[T(:,2)*100 T(:,3:4)]);
plot(ecol([1 1]),[0 1],'k','linewidth',2);
title('Transmission Probabilities for F + H_2(v=0)','fontsize',24);
set(gca,'fontsize',20);
xlabel('Collision Energy, eV');
legend('0-1(X100)','0-2','0-3','location','north');
axis([ecol([1 end]) 0 1]);

figure(2)
tricontour(mesh.p,mesh.tq,abs(psi(:,ll)).^2,25);
text(7,4.75,'FH+H','fontsize',18);
text(8,0.4,'F+H_2','fontsize',18);
title('Probability Density E_{col} = 0.197 eV','fontsize',24);
set(gca,'fontsize',20);
axis off;

figure(3);
hold on
% classically forbidden region contour
tricontour(mesh.p,mesh.tq,v,[E(ll) E(ll)]);
colormap([1,0,0])

resampledensity = 200;
% probability current interpolated to a rougher grid
[x1,y1,Jx1,Jy1] = reflux(mesh,psi(:,ll) ,m,resampledensity);
quiver(x1,y1,Jx1,Jy1);
% probability current from raw data
[jx,jy] = flux(mesh.p,mesh.t,psi(:,ll),m);

title('Probability Current E_{col} = 0.197 eV','fontsize',24);
text(7,4.75,'FH+H','fontsize',18);
text(8,0.4,'F+H_2','fontsize',18);
set(gca,'fontsize',20);
axis off;

figure(4);
curlZ   = tricurl(mesh.t,mesh.p(:,1),mesh.p(:,2),jx,jy);
tricontour(mesh.p,mesh.tq,curlZ,20);
text(7,4.75,'FH+H','fontsize',18);
text(8,0.4,'F+H_2','fontsize',18);
axis off;
title('Curl of Probability Current E_{col} = 0.197 eV','fontsize',24);