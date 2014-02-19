function mesh = fh2_mesh(n,ploton)
% set masses
mh = 1836.15264;         % mass of hydrogen in au
m = mh*[19/1.008 1 1];   % mass of [ f h h]     
vin = 0.12;              % countour of inner boundary
vout = .5;               % countour of outer boundary
rmax = 8.5;               % cutoff value for reactant 
pmax = 8;             % cutoff value for product
npar = 1.18;             % parameter to ensure correct number of points


%% potential
rhh = .3:0.05:9;
rfh = .8:0.05:9;

[x,y]=meshgrid(rfh,rhh);
pb = tojac(x(:),y(:),m);

% v=potfh2_col(x(:),y(:));
v = fh2_muck(x(:),y(:),1,[]);

v = reshape(v,size(x));
Ra = reshape(pb(:,1),size(x));
ra = reshape(pb(:,2),size(x));

%% call the mesh generator
 
mesh = meshgen(m,vin,vout,rmax,pmax,Ra,ra,v,npar,n,ploton,'fh2_');

