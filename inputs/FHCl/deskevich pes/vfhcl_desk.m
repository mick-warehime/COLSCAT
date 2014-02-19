function  [vnode,va,vc] = vfhcl_desk(uab,ubc,Av,bnd)

vnode = vfhcl_desk_col(uab,ubc);

if nargout>2
    
    % make sure the minimum of the potential is set to the minimum of the
    % reactant channel  
    vnode = vnode-min(vnode(bnd.ba));
    
    % along the boundaries
    va = vnode(bnd.ba);
    vc = vnode(bnd.bc);
end

