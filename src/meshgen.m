% generates a triangulated domain for a reactive scattering problem
function mesh = meshgen(m,vin,vout,rmax,pmax,x,y,v,spar,n,ploton,filename)


%----------------------------------------------------------------
% input parameters

% m     = [m_a m_b m_c]
% vin   = potential along inner potential boundary
% vout  = potential along oer potential boundary
% dq    = approximate side length of triangles in mesh
% rmax  = largest value of Ra (defines reactant boundary)
% pmax  = largest value of Rc (defines product boundary)
% vname = returns [uab ubc vmat]
% Av    = Coupling matrix such that v*Av = coupled potential


%% skew angle
M = mass(m);

%% get contour lines of v at vin and vout

% contour plot of potential surface
figure(100)
co = contour(x,y,v,[vout vout],'linestyle','none');
ci = contour(x,y,v,[vin vin],'linestyle','none');
close

% make sure the outer contour has the correct orientation
co = co(:,2:(co(2,1)+1))';
if co(end,2)>co(1,2);
    co = flipud(co);
end
co(end+1,:) = [rmax co(end,2)];

ci = ci(:,(ci(2,1)+3):end)';

%% fit a hyperbola to the inside and outside potential walls

mesh.pi = fminsearch(@(p) hypfit(ci,p),[2 .5 min(ci(:,1)) .2 M.alpha/2]);
mesh.po = fminsearch(@(p) hypfit(co,p),[2 .5 min(co(:,1)) .2 M.alpha./2]);

% xhyp = @(y,p) p(1)*sqrt(((y-p(4))/p(2)).^2+1)+p(3);

%%  parameters for distmesh

%fixed points in the grid
mesh.pf = hypline(mesh.pi,mesh.po,rmax,pmax,M.alpha);

% make sure pfix is not mucked up

if ~isreal(mesh.pf)
    mesh = [];
    disp('Bad mesh: check that vin is not too high')
    return
end
% bounding box
mesh.hbox = [min(co(:,1))-.5 mesh.pf(end,2); max(mesh.pf(:,1)) max(mesh.pf(:,2))+.5];

%% determine the approximate triangle sidelenght to acheive n points
area = mesharea(mesh,1e5);
mesh.dq = spar*sqrt(area/n);

%% other input data
mesh.R = [rmax,pmax]; % max value for reactants and products

%% use distmesh2d to calculate a triangulation
[p,t] = distmesh2d(@dmesh,@huniform,mesh.dq,mesh.hbox,mesh.pf,ploton,mesh);
mesh.tl = t;

%% calculate quadratic mesh
[mesh.p,mesh.t,mesh.tq] = quadmesh(p,t);

%% append boundary
mesh.bnd = meshboundary(mesh);

%% save data
fname = [filename,num2str(length(mesh.p))];
save(strcat(fname),'mesh');
close 
