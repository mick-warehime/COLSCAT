function [zx,zy] = trigradientP2(x,y,z,t)

p = [x(:) y(:)];
X = x(t)';
Y = y(t)';
Z = z(t).';

% number of triangles and points
nt = size(t,1); np = length(p);

% determine the 2nd order surface above each triangle
row = kron(1:6*nt,ones(6,1))';
col = bsxfun(@plus,kron(1:6:6*nt,ones(1,6)),(0:5)')';
coefs = sparse(row,col,[X(:).^2 X(:).*Y(:) Y(:).^2 X(:) Y(:) ones(size(X(:)))])\Z(:);
a = coefs(1:6:end);
b = coefs(2:6:end);
c = coefs(3:6:end);
d = coefs(4:6:end);
e = coefs(5:6:end);

fx = bsxfun(@times,2*a,X')+bsxfun(@times,b,Y')+kron(d,ones(1,6));
fy = bsxfun(@times,2*c,Y')+bsxfun(@times,b,X')+kron(e,ones(1,6));

% calculate triangle areas using the cross product
p21 = p(t(:,2),:)-p(t(:,1),:);
p31 = p(t(:,3),:)-p(t(:,1),:);
areas = (.5*abs(p21(:,1).*p31(:,2)-p31(:,1).*p21(:,2)))';

% map matrix from triangles to nodes
M = sparse(repmat((1:nt)',6,1),t(:),1,nt,np);
 
% weight the partial derivative (a or b coefficient) by the area of each
% triangle and sum for a given node. divide by the total area 
fxa = bsxfun(@times,fx,areas');
fya = bsxfun(@times,fy,areas');

zx =  full(diag(sparse(t,t,fxa,np,np)))./(areas*M)';
zy =  full(diag(sparse(t,t,fya,np,np)))./(areas*M)';

