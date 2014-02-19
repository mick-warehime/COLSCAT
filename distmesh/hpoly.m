function h = hpoly(p,mesh,vfun)

% inputs
% p = mesh points from distmesh2d
% pv = pfixed from input to distmesh2d
% f = TriScatteredInterp(uab,ubc,v) from ab initio data 
%     (f must be in mass scaled jacobi coordinates)

% evaluate the potential at the mesh points
v  = vfun(p(:,1),p(:,2));

% use the potential values to determine the 'h' function for distmesh2d
% larger potential value --> larger mesh in that region
% smaller potential value --> smaller mesh in that region
h = 2-exp(-v/max(v)).^3;

