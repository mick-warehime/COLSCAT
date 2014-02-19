% resamples the flux inside the boundary. 

% input

% p  - original points in the mesh
% pb - polygonal domain to the boundary
% pf - fixed points on the boundary used in distmesh2d
% jx - flux in the x dir
% jy - flux in the y dir
% dx - new x dir step size
% dy - new y dir step size

% output

% jxn - flux in x dir sampled at new points
% jyn - flux in y dir sampled at new points

function [xn yn jxn jyn] = reflux(p,pb,jx,jy,dq)

% tol remove points that are less that tol% of the max(abs(jx^2+jy^2))
tol = 1e-5;

% create the sampling functions
JX = TriScatteredInterp(p(:,1),p(:,2),jx);
JY = TriScatteredInterp(p(:,1),p(:,2),jy);

% use mesh2d to calculate new smaller mesh
warning off
hdata.fun = @hfun;
hdata.args = {dq};
options.output= false;
[pn,t]=mesh2d(pb,[],hdata,options);

% store points
xn = pn(:,1);
yn = pn(:,2);

% sample at new points
jxn = JX(xn,yn);
jyn = JY(xn,yn);

% indices to remove values that are nearly zero
J = jxn.^2+jyn.^2;
notok = J<(max(J)*tol);

% delete points
jxn(notok) = [];
jyn(notok) = [];
xn(notok) = [];
yn(notok) = [];

end % reflux()

function h = hfun(x,y,dq)


h = dq*ones(size(x));

end % hfun

