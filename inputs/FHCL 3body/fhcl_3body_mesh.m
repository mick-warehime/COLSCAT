function mesh = fhcl_3body_mesh(Av,n,vinput,ploton)

% set masses
mh = 1836.15264;                    % mass of hydrogen in au
m  = mh*[19.00 1.008 35.453]/1.008; % mass of [F,H,Cl] in amu
vin = 0.06;                         % countour of inner boundary
vout = 0.3;                         % countour of outer boundary
rmax = 13;                          % cutoff value for reactant
pmax = 13;                          % cutoff value for product
npar = 1.7;                         % larger n = fewer points
 
%% potential
[x,y] = meshgrid(.6:.05:10,1:.05:12);
pb = tojac(x(:),y(:),m);

v =  feval(vinput,x(:),y(:),Av,[]);

vj = reshape(v(:,1),size(x));
xj = reshape(pb(:,1),size(x));
yj = reshape(pb(:,2),size(x));

%% call the mesh generator
mesh = meshgen(m,vin,vout,rmax,pmax,xj,yj,vj,npar,n,ploton,'fhcl_');

