% resample the flux with fewer points

function [x_new,y_new,Jx_new,Jy_new] = reflux(mesh,psi,m,n)

[Jx,Jy] = flux(mesh.p,mesh.t,psi,m);

% range of points in mesh
box = [range(mesh.p(:,1))' range(mesh.p(:,2))'];

% area of mesh
marea = mesharea(mesh,1e5);

% approximate number of points in the box to get the desired number in the mesh
dn = round(sqrt(n*marea));
[x,y] = meshgrid(box(1,1):diff(box(:,1))/(dn-1):box(2,1),box(1,2):diff(box(:,2))/(dn-1):box(2,2));
x(1:2:end,:) = x(1:2:end,:)+ diff(box(:,1))/(2*(dn-1)); % shift every other row
in = inpolygon(x(:),y(:),mesh.bnd(:,1),mesh.bnd(:,2));
x_new = x(in);
y_new = y(in);

% interpolate to new points
Fx = scatteredInterpolant(mesh.p(:,1),mesh.p(:,2),Jx);
Fy = scatteredInterpolant(mesh.p(:,1),mesh.p(:,2),Jy);
Jx_new = Fx(x_new,y_new);
Jy_new = Fy(x_new,y_new);

% only return values of the flux that are not neglible
in2 = sqrt(Jx_new.^2+Jy_new.^2)>1e-3*max(sqrt(Jx_new.^2+Jy_new.^2));
x_new = x_new(in2);
Jx_new = Jx_new(in2);
y_new = y_new(in2);
Jy_new = Jy_new(in2);