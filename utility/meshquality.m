% returns the 'quality' for each element in a triangular mesh
% from persson/strang paper point 7: mesh quality:

% q = (b+c-a)*(c+a-b)*(a+b-c)/(a*b*c)
% equilateral triangles have q=1 one
% degenerate triangles (zero area) have q=0
% if q>0.5 for all elements, the mesh should be good

function [q,r] = meshquality(p,t)

% calculate the side lengths
a = sqrt( sum( (p(t(:,1),:) - p(t(:,2),:) ) .^2 ,2) ) ;
b = sqrt( sum( (p(t(:,2),:) - p(t(:,3),:) ) .^2 ,2) ) ;
c = sqrt( sum( (p(t(:,3),:) - p(t(:,1),:) ) .^2 ,2) );

% return the mesh quality
q = (b+c-a).*(c+a-b).*(a+b-c)./(a.*b.*c);

% circumscribed radius
r = a.*b.*c./sqrt((a+b+c).*(b+c-a).*(c+a-b).*(a+b-c));
