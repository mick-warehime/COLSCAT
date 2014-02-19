% return the divergence of the input wave function

function div = divpsi(p,t,m,psi)

if ~iscell(psi)
    psi = {psi};
end

div = cell(length(psi),1);

for j=1:length(psi)
    
    for k = 1:size(psi{1},2)
        [Jx,Jy] = flux(p,t,psi{j}(:,k),m);
        div{j}(:,k) = tridivergence(p,t,[Jx,Jy]);
    end
end