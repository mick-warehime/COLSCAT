% resample the flux on a given triangulation

function [x y JX_coarse JY_coarse] = coarsenflux(mesh,JX,JY,dx)

p = mesh.p;

% create a new coarer mesh grid
[xm ym] = meshgrid(min(p(:,1)):dx:max(p(:,1)),min(p(:,2)):dx:max(p(:,2)));

% remove points outside the boundary of the mesh
d = dmesh([xm(:) ym(:)],mesh.pf,mesh.pi,mesh.po);
x = xm(d<0); y = ym(d<0);

% interpolate to new points
JX_coarse = cell(length(JX),1);
JY_coarse = cell(length(JX),1);
for j=1:length(JX)
    Fx1 = TriScatteredInterp(p,JX{j}(:,1));
    Fx2 = TriScatteredInterp(p,JX{j}(:,2));
    Fy1 = TriScatteredInterp(p,JY{j}(:,1));
    Fy2 = TriScatteredInterp(p,JY{j}(:,2));
    
    JX_coarse{j} = [Fx1([x,y]) Fx2([x y])];
    JY_coarse{j} = [Fy1([x,y]) Fy2([x y])];
    
end

