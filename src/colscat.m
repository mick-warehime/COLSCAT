% Quantum Coupled-state Scattering
function [R,T,psi,mesh,jac,viba,vibc,vnode,bnd]  = colscat(m,E,initState,serf,mesh,Av,vinput)

%% number of electronic surfaces
if size(Av,2) == 1; ns = 1; else ns = 2;end

%% get a mass structure
M = mass(m);

%% get the indices of the boundaries
[mesh,np,bnd] = boundary(mesh);

%% Boundary Jacobi Coordinates
jac = jacobi(mesh.p,bnd,M);

%% spline potential at nodes and the center of each triangle
pb = tobond(mesh.p,m);
[vnode,va,vc] = feval(vinput,pb(:,1),pb(:,2),Av,bnd);

%% Calculate vibrational wave functions
viba = cell(ns,1); vibc = cell(ns,1);
for j=1:ns
    viba{j} = vibfemP2(jac.ra,va(:,j),M.mubc);
    vibc{j} = vibfemP2(jac.rc,vc(:,j),M.muab);
end

%% FEM MATRICES (2nd order 2D)
[T,V,O] = femP2(mesh,vnode,ns);

%% Scattering probability as a function of energy
le = length(E); S = cell(le,1); PSI = cell(le,1);
parfor jj = 1:le
    [PSI{jj},S{jj}] = solver(E(jj),T,V,O,bnd,M,viba,vibc,initState,serf,ns,np);
end

% post process the parallel data
[R,T,psi] = stort(S,PSI,le,ns,np);

