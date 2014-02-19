% generate a mesh for the HH2 reaction
function mesh = hh2_mesh(n,ploton)
% set masses
mh = 1836.15264;       % mass of hydrogen in au
m = [mh mh mh];        % 
vin = 0.1;             % countour of inner boundary
vout = 0.6;            % countour of outer boundary
rmax = 6;              % cutoff value for reactant 
pmax = 6;              % cutoff value for product
npar = 1.5;            % parameter to ensure correct number of points

%% potential
vec = 0.25:0.05:8.9;

[x,y]=meshgrid(vec,vec);
pb = tojac(x(:),y(:),m);

v=vhh2(x(:),y(:));

v = reshape(v,size(x));
Ra = reshape(pb(:,1),size(x));
ra = reshape(pb(:,2),size(x));

%% call the mesh generator
mesh = meshgen(m,vin,vout,rmax,pmax,Ra,ra,v,npar,n,ploton,'hh2_');

