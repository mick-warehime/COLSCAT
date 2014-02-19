function [T,V,O] = femP2(mesh,v,ns)
p = mesh.p;         % nodal points
t = mesh.t(:,1:3); % linear mesh
tq = mesh.t;        % full quadratic mesh

% compute the area of each triangle (half of the crossproduct)
r1 = p(t(:,2),:)-p(t(:,1),:);
r2 = p(t(:,3),:)-p(t(:,1),:);
areas = .5*abs(r1(:,1).*r2(:,2)-r2(:,1).*r1(:,2));

%% Code to calculate A and D adapted from Jan Valdman's File Exchange Code
% "Fast assembly of stiffness and matrices in finite element method using nodal elements"

Xscalar=kron(ones(1,6),tq); 
Yscalar=kron(tq,ones(1,6));

% get sub matrices
[k0,k1,k2] = p2matrix(p,t);
 
% overlap matrix: int \phi_i \phi_j
O=sparse(Xscalar,Yscalar,bsxfun(@times,k2,2*areas));

% kinetic energy matrix: int \del \phi_i \del \phi_j
T=sparse(Xscalar,Yscalar,bsxfun(@rdivide,k0,2*areas));

if ns == 2   
    T = blkdiag(T,T);
    O = blkdiag(O,O);    
end

%% Calculate the potential energy matrix element integrals

% to integrate (phi_i*phi_j*v) over the standard triangle we can simply
% take the correct weighted sum of the potential points
V = cell(ns);
for i=1:ns
    for j=1:ns
        vv = v(:,i + j -1);
        V{i,j} = sparse(Xscalar,Yscalar,bsxfun(@times,vv(tq)*k1,areas)); % instead of 2*areas the factor of two is used to simplify c matrix in p2matrix
    end
end
V = cell2mat(V);
