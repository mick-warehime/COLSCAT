void
%% Input

% set masses
me = 1822.88862599;                 % electron mass
m = [18.998, 1.008, 35.453]*me;     % mass of [F,H,Cl] in amu
E = 0.008:.0001:.06;               % energy range hartree
initState = 0;                          % initial state
serf = 1;                           % initial surface
  
%% generate mesh or load mesh

% create new mesh
% n=2500;
% mesh = fhcl_mesh(n,1);

% or load a previously calculated mesh
load fhcl_desk_11705

 

%% Call main Szcattering Routine
[R,T,psi,mesh,jac,viba,vibc,vnodes,bnd] = colscat(m,E,initState,serf,mesh,1,'vfhcl_desk');
 
% collision energy in eV
ecol = (E-viba{1}.e(1))*27.211;

% plot reflection coefficients for first three vibrational states
figure(1)
plot(ecol,R(:,1:5))
legend('0-0','0-1','0-2')
title('Reflection Probabilities for F+ HCl(v=0)','fontsize',24)
set(gca,'fontsize',20)
axis([ecol([1 end]) 0 1])
set(gca,'xtick',0:.25:1.5)
xlabel('Collision Energy, eV')

% plot reflection coefficients for first three vibrational states
figure(2)
plot(ecol,T(:,1:5))
legend('0-0','0-1','0-2')
title('Transmission Probabilities for F + HCl(v=0)','fontsize',24)
set(gca,'fontsize',20)
axis([ecol([1 end]) 0 1])
set(gca,'xtick',0:.25:1.5)
xlabel('Collision Energy, eV')