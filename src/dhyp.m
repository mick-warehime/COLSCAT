% given the parameters of a rotated hyperbola return something awesome

function d = dhyp(p,par)
% model the x values using a hyperbolic fit ((x-xc)/a)^2-(y-yc)^2/b^2=1
% par = [a b xc yc]

% rotate the points
pr = protate(p,par(:,5));    

% the implicit function of a hyperbola is used to calculate distance
d = 1- ( ((pr(:,1)-par(3))/par(1)).^2-((pr(:,2)-par(4))/par(2)).^2);

xin = pr(:,1)<par(3);

d(xin) = max([d(xin) (p(xin,1)-par(3)).^4],[],2);
% % sometimes you just want the xvalues for a given set of y values
% x =  par(1)*sqrt(1+(pr(:,2)-par(4)).^2/par(2)^2)+par(3);
% 
% xy = protate([x,pr(:,2)],-a);
