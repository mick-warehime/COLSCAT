function [PSI,S] = solver(E,T,V,O,B,M,viba,vibc,initState,serf,ns,np)

% account for matlabs vector index v=[0,1,2,3] --> state=[1,2,3,4]
initState = initState+1;

% define the G matrix
A = T+2*M.mu*(V-E*O);

% extend the system using boundary conditions
[Q,b] = extendsys(A,viba,vibc,initState,serf,M,B,E,np,ns);

% simple preconditioner
% C = diag(max(abs(Q.').^2));

% solve linear system
U = Q\b;

% wave function
psi = full(U(1:(np*ns),:));

% scattering coefficients
s = abs(full(U((ns*np+1):end,:))).^2;

%% extract and sort scattering coefficients
nba = length(B.ba); nbc = length(B.bc);
S = cell(ns,1); PSI = zeros(np,ns);
for jj = 1:ns
    % indices for a and c channels on given surface
    ain = (jj-1)*(nba+nbc);
    cin = nba+(jj-1)*(nba+nbc);  

    % store probabilities in each channel on each surface
    S{jj}.a = s((1:nba)+ain,:);
    S{jj}.c = s((1:nbc)+cin,:)*(M.mua/M.muc);
    PSI(:,jj) = psi((1:np)+(jj-1)*np,:);
end

