% return the distance to the reaction path

function h = hmesh(p,mesh)

mh = 1836.15264;       % mass of hydrogen in au
m  = [mh,mh,mh]; 
Av = 1;

% convert to bond coordinates
msj = tobond(p,m);

% distance to the reaction path
v = vhh2(msj(:,1),msj(:,2),Av); 

% exponentiate
h = exp(v-min(v)+1);
