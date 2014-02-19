% calculate the flux of the wavefunction defined pointwise over a given
% triangulation

function [Jx,Jy] = flux(p,t,psi,m)

% preallocate
Jx = zeros(size(psi));
Jy = zeros(size(psi));

% reduced mass of the triatomic system in mass scaled jacobi coordinates
M = mass(m);

% determine if the triangulation is linear or quadratic
[~,trisize] = size(t);
% get the partial derivatives of the wave function
if trisize==3    
    [Fx,Fy] = trigradientP1(p(:,1),p(:,2),psi,t);
 elseif trisize==6
    [Fx,Fy] = trigradientP2(p(:,1),p(:,2),psi,t);
end

for j=1:size(psi,2)  
        
    % calculate the flux
    Jx(:,j) = -(1i/(2*M.mua))*((psi(:,j)').'.*Fx(:,j)-psi(:,j).*(Fx(:,j)').');
    Jy(:,j) = -(1i/(2*M.mua))*((psi(:,j)').'.*Fy(:,j)-psi(:,j).*(Fy(:,j)').');
    
end