% matrices used to extend the scattering problem
function [B,F,I,bo,fo] = extendmat(viba,vibc,initState,serf,M,bnd,E,np,ns,currsrf)

nstate = length(initState);
nba = length(bnd.ba); 
nbc = length(bnd.bc);

% get the wave vectors
ka = sqrt(2*M.mua*(E-viba{currsrf}.e')); na = length(ka);
kc = sqrt(2*M.muc*(E-vibc{currsrf}.e')); nc = length(kc);
ki = zeros(1,nstate);
for jj=1:nstate
    ki(jj) = sqrt(2*M.mua*(E-viba{serf}.e(initState(jj))'));
end 

% multiply the integrated wavefunctions by the appropriate derivatives
Ba = bsxfun(@times,viba{currsrf}.iwf,1i*sqrt(ka))/M.lama^2;
Bc = bsxfun(@times,vibc{currsrf}.iwf,1i*sqrt(kc))/M.lamc^2;

B = sparse(repmat(bnd.ba,1,na),repmat(1:na,nba,1),Ba,np,na+nc)+...
    + sparse(repmat(bnd.bc,1,nc),repmat(1:nc,nbc,1)+na,Bc,np,na+nc);

% I is used to contruct the boundary term
I = sparse(1:(nba+nbc),[bnd.ba;bnd.bc],ones(nba+nbc,1),nba+nbc,np);

% the H matrix (wave function boundary conditions)
Fa = bsxfun(@times,viba{currsrf}.wf,sqrt(1./ka));
Fc = bsxfun(@times,vibc{currsrf}.wf,sqrt(1./kc));
F = blkdiag(Fa,Fc);

if serf == currsrf    
    % the b matrix (wave function boundary conditions)        
    fio =  bsxfun(@times,viba{currsrf}.wf(:,initState),sqrt(1./ki));
    
    row_fo = repmat(((1:nba)+(currsrf-1)*(nba+nbc))',1,nstate);
    col_fo = kron(ones(nba,1),1:nstate);
    fo  = sparse(row_fo,col_fo,fio,(nba+nbc)*ns,nstate);

    % multiply the integrated wavefunctions by the appropriate derivatives
    bio = bsxfun(@times,viba{currsrf}.iwf(:,initState),-1i*sqrt(ki))/M.lama^2;
    
    row_bo = repmat((bnd.ba+np*(currsrf-1)'),1,nstate);
    col_bo = kron(ones(nba,1),1:nstate);
    bo  = sparse(row_bo,col_bo,bio,ns*np,nstate);    
else
    bo = sparse(ns*np,nstate);
    fo = sparse((nba+nbc),nstate);
end
 