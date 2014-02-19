% interpolate data of the Muckerman 5 LEPS potential

function [vnode,va,vc] = fh2_muck(x,y,Av,bnd)

% load muckerman data
m5data = load('muckfh2');

% meshgrid of data points
[rfh,rhh] = meshgrid(m5data.ua,m5data.uc);

% interpolate data to mesh points
vnode = interp2(rfh,rhh,m5data.v,x,y,'cubic');

if nargout>1
    

    % make sure the minimum of the potential is set to the minimum of the
    % reactant channel  
    vnode = vnode-min(vnode(bnd.ba));
    
    % along the boundaries
    va = vnode(bnd.ba);
    vc = vnode(bnd.bc);
    
end