% reorganize the scattering coefficients

function [R,T,psi] = stort(S,PSI,le,ns,np)

% number of output states
n = 10;


% sort the probabilities and wave functions for single states
if ns == 1
    R = zeros(le,n);
    T = zeros(le,n);
    psi = zeros(np,le);
    for jj = 1:le
        for kk = 1:n
            R(jj,kk) = S{jj}{1}.a(kk);
            T(jj,kk) = S{jj}{1}.c(kk);
        end
        psi(:,jj) = PSI{jj};
    end
    
elseif ns == 2
    
    R = cell(ns,1);
    T = cell(ns,1);
    psi = cell(ns,1);
    for jj=1:ns
        for kk = 1:le
            for ll = 1:n
                R{jj}(kk,ll) = S{kk}{jj}.a(ll);
                T{jj}(kk,ll) = S{kk}{jj}.c(ll);
            end
            psi{jj}(:,kk) = PSI{kk}(:,jj);            
        end
        
    end
    
    
end

