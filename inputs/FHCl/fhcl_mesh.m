function mesh = fhcl_mesh_desk(n,ploton)

% set masses
me = 1822.88862599;
m = [18.998, 1.008, 35.453]*me; % mass of [F,H,Cl] in au

vin = 0.06;             % countour of inner boundary
vout = 0.3;             % countour of outer boundary
rmax = 13;              % cutoff value for reactant
pmax = 13;              % cutoff value for product
npar = 1.7;             % larger n = fewer points
 
%% potential

[x,y] = meshgrid(1:.01:5,2:.01:5);
pb = tojac(x(:),y(:),m);

vj = feval('vfhcl_desk',x,y,1,[]);
xj = reshape(pb(:,1),size(x));
yj = reshape(pb(:,2),size(x));

%% call the mesh generator
mesh = meshgen(m,vin,vout,rmax,pmax,xj,yj,vj,npar,n,ploton,'fhcl_desk_');
 
