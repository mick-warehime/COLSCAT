% convert from bond coordinates to mass scaled jacobi coordinates

function p = tojac(x,y,m)

M = mass(m);

% convert to msjacobi coordinates
p(:,1) = (x+y*M.mc/(M.mb+M.mc))*M.lama;
p(:,2) = y/M.lama;