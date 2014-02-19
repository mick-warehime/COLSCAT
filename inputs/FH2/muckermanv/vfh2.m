function [vnode,va,vc] = vfh2(x,y,Av,bnd)


vnode = potfh2_col(x(:),y(:));


if nargout>1
    

    % make sure the minimum of the potential is set to the minimum of the
    % reactant channel  
    vnode = vnode-min(vnode(bnd.ba));
    
    % along the boundaries
    va = vnode(bnd.ba);
    vc = vnode(bnd.bc);
    
end