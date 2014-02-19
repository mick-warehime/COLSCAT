function [vnode,va,vc] = vhh2(x,y,Av,B)

% load mielke potential data
vdata = importdata('hh2_data.txt');

ncol = length(unique(vdata(:,1)));
nrow = length(unique(vdata(:,2)));

uab = reshape(vdata(:,1),nrow,ncol);
ubc = reshape(vdata(:,2),nrow,ncol);
vv = reshape(vdata(:,3),nrow,ncol);

vnode = interp2(uab,ubc,vv,x,y,'cubic');


if nargout>1
    
    % zero of the potential is set to the minimum of the reactant channel  
    vnode = vnode-min(vnode(B.ba));
    
    % potential along the boundaries to calculate vibrational bound states
    va = vnode(B.ba);
    vc = vnode(B.bc);
end