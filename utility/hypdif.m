% return the values of the implicit hyperbolic function and its first and
% second derivatives

function [f fx fxx fy fyy] = hypdif(p,par)

% model the x values using a hyperbolic fit ((x-xc)/a)^2-(y-yc)^2/b^2=1
% par = [a b xc yc]
pr = protate(p,par(:,5));    

f = ((pr(:,1)-par(3))/par(1)).^2-((pr(:,2)-par(4))/par(2)).^2;

fx = 2*(pr(:,1)-par(3))/par(1)^2;

fy = -2*(pr(:,2)-par(4))/par(2)^2;

fxx = 2/par(1)^2;

fyy = -2/par(2)^2;