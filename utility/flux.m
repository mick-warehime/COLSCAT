% calculate the flux of the wavefunction defined pointwise over a given
% triangulation

function [Jx,Jy] = flux(p,t,psi,m)

% preallocate
Jx = zeros(size(psi));
Jy = zeros(size(psi));

% reduced mass of the triatomic system in mass scaled jacobi coordinates
mu = sqrt(m(1)*m(2)*m(3)/sum(m));

% get the partial derivatives of the wave function
[Fx,Fy] = trigradient(p,t,psi);

for j=1:size(psi,2)  
        
    % calculate the flux
    Jx(:,j) = (psi(:,j)').'.*Fx(:,j)-psi(:,j).*(Fx(:,j)').';
    Jy(:,j) = (psi(:,j)').'.*Fy(:,j)-psi(:,j).*(Fy(:,j)').';
    
end

Jx = -(1i/(2*mu))*Jx;
Jy = -(1i/(2*mu))*Jy;