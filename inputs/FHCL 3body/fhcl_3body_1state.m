void
%% Input

% set masses
me = 1822.88862599;                 % electron mass
m = [18.998, 1.008, 35.453]*me;     % mass of [F,H,Cl] in amu
E = 0.0067:.0001:.04;               % energy range hartree
state = 0;                          % initial state
serf = 1;                           % initial surface

vinput = 'vfhcl_3body';
 Av = [2/3;1/3;  -1]; 
Av = [1; 0; 0];

%% generate mesh or load mesh
% load mesh
n=2500;
mesh = fhcl_3body_mesh(Av(:,1),n,vinput,1);

% load fhcl_desk_30802
%% Call main Szcattering Routine
tic
[R,T,psi,mesh,jac,viba,vibc,v,bnd] = colscat(m,E,state,serf,mesh,Av,vinput);
toc

ecol = (E-viba{1}.e(1))*27.211;

figure(1)
plot(ecol,T)


 
