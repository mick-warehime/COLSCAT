% looks through the generated mesh points for the points along the product
% and react boundaries

function [mesh,np,B] = boundary(mesh)

% mesh.pf(1,1:2)
%\
% \      product boundary (c channel)
%  \
%   \mesh.pf(2,1:2)
%   /r
%  /
% /
% ------mesh.pf(3,1:2)
%      |
%      |     reactant boundary (a channel)
%      |
% ------mesh.pf(4,1:2)

% ------------------------------------------------------------------------
% tolerance (used to check how far away from the boundary a given point is)
tol = 1e-5;

% fixed points on boundary
pf = mesh.pf;

p = mesh.p;
% reactant boundary
% all points with x values within tol of pf(3,1) and y values less than pf(3,2)
ba = find( abs(p(:,1)-pf(3,1))<tol & p(:,2)<=(pf(3,2)+tol) );
p(ba,1) = mean(p(ba,1));

% the product boundary goes through the pf points 1 and 2
% determine the slope and intercept of this line
m = (pf(1,2)-pf(2,2))/(pf(1,1)-pf(2,1));
b = pf(1,2)-m*pf(1,1);

% find all points within tol of the line between pf(1,:) and pf(2,:)
qline = p(:,1)*m+b;
bc = find((abs(p(:,2)-qline)-tol)<0 & p(:,2) >= pf(3,2));

% fix product points to a single line
p(bc,2) = m*p(bc,1)+b;

% sort the boundaries in order of increasing ra/rc
[~,in1] = sort(p(ba,2));
B.ba = ba(in1);
[~,in2] = sort(p(bc,1));
B.bc = bc(in2);

% number of points
np = length(p);

% store the mesh with shifted boundaries
mesh.p = p;