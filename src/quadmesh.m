% lin2quad takes a linear delauny triangulation as input (3 nodes per triangle)
% and returns a quadratic triangulation (6 nodes per triangle)


function [pq,tq,tr] = quadmesh(p,t)

[edges ,t2e] = tri2edge(t);
tq = [t t2e+length(p)];
pq = [p; (p(edges(:,1),:)+p(edges(:,2),:))/2];
tr = [tq(:,[2 4 5]); tq(:,[6 5 3]); tq(:,[4 5 6]); tq(:,[1 4 6])];