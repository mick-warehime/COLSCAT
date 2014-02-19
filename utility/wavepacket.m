% wave packet

function [F Ce ]= wavepacket(mesh,m,psi,E,t)
psi = cell2mat(psi');
p=mesh.p;
M = mass(m);

xo = 13.5; sx = .5;
yo = 1.275; sy = .25;

wp = exp(-((p(:,1)-xo)/sx).^2-((p(:,2)-yo)/sy).^2).*exp(1i*sqrt(2*M.mua*0.02)*p(:,1));

Ce = psi\wp;
tmat = exp(-1i*bsxfun(@times,t,E'));

psit=psi*bsxfun(@times,Ce,tmat);
% psit = bsxfun(@rdivide,psit,max(psit));
figure('Renderer','zbuffer')
F(length(t)) = struct('cdata',[],'colormap',[]);
set(gca,'NextPlot','replaceChildren');

for i = 1:length(t)
    
   trisurf(mesh.t,p(:,1),p(:,2),psit(:,i));
%    axis([min(p(:,1)) max(p(:,1)) min(p(:,2)) max(p(:,2)) -1 1])
   drawnow
    
    
end

