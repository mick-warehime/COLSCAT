% Fill in the submatrices of EQN. 20 by extending the linear system
function [Q,b] = extendsys(A,viba,vibc,initState,serf,M,bnd,E,np,ns)

% single state
if ns == 1
    [B,F,I,bo,fo] = extendmat(viba,vibc,initState,1,M,bnd,E,np,ns,1);
    Q = [A -B; I,-F];
    b = [bo; fo];   
     
% two states
elseif ns == 2
    
    [B1,F1,I1,bo1,fo1] = extendmat(viba,vibc,initState,serf,M,bnd,E,np,ns,1);
    [B2,F2,I2,bo2,fo2] = extendmat(viba,vibc,initState,serf,M,bnd,E,np,ns,2);
    
    Q = [   A          -blkdiag(B1,B2);
        blkdiag(I1,I2) -blkdiag(F1,F2)];

    % check to see if initial state boundary vector should be on surface
    % one or surface two
    if serf==1; b = [bo1; fo1]; else b = [bo2; fo2]; end
    
end




