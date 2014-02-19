void
%% Input

% set masses
me = 1822.88862599;                 % electron mass
m = [18.998, 2*1.008, 35.453]*me;     % mass of [F,H,Cl] in amu
E = 0.0067:.00001:.03;               % energy range hartree
coupling = 'full_j';                % chose scattering state(s), sigma, pi or full
state = 0;                          % initial state
serf = 1;                           % initial surface

vinput = 'vfhcl_3body';
 
  Av = [2/3,  sqrt(2)/3, 2/3;   % vsig
        1/3,  -sqrt(2)/3,1/3;                         % vpi
        -1,   0,         2];                          % vso
 

%% generate mesh or load mesh
% load mesh
n=17500;
mesh = fdcl_mesh(Av(:,1),n,vinput,1);
%  
% load fdcl_50552
%% Call main Szcattering Routine
tic
[R,T,psi,mesh,jac,viba,vibc,v,B] = colscat(m,E,state,serf,mesh,Av,vinput);
toc

figure(1)
plot(E,R{1})
figure(2)
plot(E,R{2})
figure(3)
plot(T{1})
figure(4)
plot(E,T{2})
 