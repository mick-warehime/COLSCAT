%% Vibrational Wavefunctions using P2FEM in 1D 

% INPUTS:
%
% r - vector of position values
% v - vector of potential as a function of r
% m - mass or reduced mass of sysmte

% OUTPUTS
%
% vib.wf - wave functions chi_n(r)
% vib.e  - eigen energies e_n
% vib.iwf  - integrated wave functions \int_\r \phi_i(r) \chi_n(r)


%% numerical approximation

function vib = vibfemP2(r,v,mu)

% check to make sure r and v are column vectors
if size(r,1) == 1; r = r';end
if size(v,1) == 1; v = v';end

% r vector length and step size
lnr = length(r);
dr = diff(r);
dr = dr(1:2:end);

% sparse index labels
in = bsxfun(@plus,[1 2 3],(0:2:(lnr-3))');
xin = kron([1 1 1],in);
yin = kron(in,[1 1 1]);

% weight factors (SEE P2integrals.m)
k0 =[ 7, -8, 1, -8, 16, -8, 1, -8, 7]/6;
k1 = [ 39, 20, -3, 20,  16, -8, -3, -8, -3;...
    20, 16, -8, 16, 192, 16, -8, 16, 20;...
    -3, -8, -3, -8,  16, 20, -3, 20, 39]/210;
k2 =[ 4, 2, -1, 2, 16, 2, -1, 2, 4]/15;

%% del phi matrix

K0 = sparse(xin,yin,bsxfun(@rdivide,k0,dr));

%% potential matrix

K1 = sparse(xin,yin,bsxfun(@times,v(in)*k1,dr));

%% overlap matrix

K2 = sparse(xin,yin,bsxfun(@times,k2,dr));

%% Calculate the energies and wave functions

% solve the Ax = 0 FEM problem to determine lnr wavefuncs and energies
[wf,e] = eigs(blkdiag(K0/(2*mu)+K1,1),blkdiag(K2,1),lnr,'sm');

% make sure lowest energy state is first
% (eigs(A,B,numstates,'sm') returns smallest energy in last column
[EWF,eindex]  = sort(diag(e));
wf = wf(1:(lnr),eindex);

%% Integrate the vibrational wavefuncitons along the boundary
% Fij = int(phi_i chi_j)
% expand Chi_j = \sum_i c_i \phi_i
% then Fi = D*wf_i

% F = K2*wf;

%% store results

vib.r = r;
vib.wf = bsxfun(@rdivide,wf,sqrt(trapz(r,wf.^2))); % normalize wavefunctions
vib.e = EWF;
vib.iwf = K2*wf;

