function [vnode,va,vc] = vfhcl_3body(x,y,Av,bnd)

% sigma and pi potentials using 3body fit
[mhf,mhcl,D] = importpotential('fhcl_parameters_3body_new.txt');

% v3 = [V_sig V_pi]
v3 = v3body(x(:),y(:),mhf,mhcl,D);

% spin orbit constant
vso = fhcl_so_approx(x(:),y(:),100);

% [vsig vpi vso] * Av ===> Av selects the basis
vnode = [v3 vso]*Av;

if nargout>1
    
    % along the boundaries
    diagstates = 1;
    if size(Av,2) >1
        diagstates(2) = 3;
    end
    
    % make sure the minimum of the potential is set to the minimum of the
    % reactant channel  
    vnode(:,diagstates) = vnode(:,diagstates)-min(vnode(bnd.ba,1));
    
    % along the boundaries
    va = vnode(bnd.ba,diagstates);
    vc = vnode(bnd.bc,diagstates);
end