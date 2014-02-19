% distance to the reactant and product boundary

function d = dmesh(p,mesh)

%% distance to the hyperbolic boundaries
dhi = dhyp(p,mesh.pi);
dho = dhyp(p,mesh.po);

%% distance to the product/reactant boundaries

% indices for three sections of points
rin = p(:,2)<=mesh.pf(3,2);
pin = p(:,2)>=mesh.pf(2,2);
min = ~rin & ~pin; 

% distance to reactant boundary
dr  = p(:,1)-mesh.pf(3,1);

% distance to middle section (avoids reactant boundary affecting product)
% formulated to check if x values are to the left of the boundary line
mm = (mesh.pf(3,1)-mesh.pf(2,1))/(mesh.pf(3,2)-mesh.pf(2,2));
bm = mesh.pf(3,1)-mm*mesh.pf(3,2);
dm = p(:,1)-(mm*p(:,2)+bm);

% distance to product boundary
mp = (mesh.pf(1,2)-mesh.pf(2,2))/(mesh.pf(1,1)-mesh.pf(2,1));
bp = mesh.pf(1,2)-mp*mesh.pf(1,1);
dp = p(:,2)-(mp*p(:,1)+bp);

% use patchwork distance function for boundary section
drmp(rin) = dr(rin);
drmp(min) = dm(min);
drmp(pin) = dp(pin);

%% combine all distance functions
d = max([-dhi max([ dho drmp' ],[],2)],[],2);

%% diagnostic plots

% figure(1)
% in = d<0;
% plot(p(in,1),p(in,2),'ko',mesh.pf(:,1),mesh.pf(:,2),'ro-')
% 
% figure(2)
% ini = dhi<0;
% ino = dho<0;
% plot(p(ino,1),p(ino,2),'ro',p(ini,1),p(ini,2),'ko')
% 
% figure(3)
% inrmp = drmp<0;
% plot(p(inrmp,1),p(inrmp,2),'ro')


