% calculate all the integrals for the 1d fem calculation used to determine
% the vibrational wave functions and energies

%  [(-1/(2mu))d^2/dx^2  + v(x)]\psi = E \psi 
%  
%  premultiply by xi and integrate over boundry 
%  use int. by parts for the integral
% \int \xi d^2/dx^2 \psi = \xi \psi  - \int d/dx \xi d/dx \psi
%
%  substituting with hat functions the problem is written
%%
% $ \frac{1}{2\mu}\int \frac{d}{dx} \phi_i \frac{d}{dx} \phi_j + \int \phi_i v \phi_j = E \int \phi_i \phi_j $

%%
%   in matrix notation
%
%   [1/(2mu)*D + V]c = E*O
%%
% $x^2+e^{\pi i}$

%     hat functions and their overlap with neighbors only 3 nonzero elements
%     
%    /|\  /|\  /|\
%   / | \/ | \/ | \ 
%  /  | /\ | /\ |  \
% /   |/  \|/  \|   \
%     o    o    o
%     i    j    k
%   
%         x=0 

% Each hat function is defined as follows:
%
% $$\phi_j = \left{  \begin{array}{c} 1+\frac{x}{h_{ij}} \\ 1-\frac{x}{h_{jk}} \end{array}   \right.$$
%
%    where $h_{lm}$ = x_m - x_l
% 

%% define basis functions
void

% set x_j = 0

syms hij hjk x vi vj vk 

phi_i = -x/hij;

phi_jl = 1+x/hij;
phi_jr = 1-x/hjk;

phi_k = x/hjk;


%%  overlap of the derivatives

% D_ij = \int d/dx \phi_{i} d/dx \phi_{j}

% off diagonal terms

D_ij = int( diff(phi_i,x)*diff(phi_jl) ,x ,-hij ,0);    % = -1/hij

D_jk = int( diff(phi_k,x)*diff(phi_jr) ,x ,0 ,hjk);     % = -1/hjk

% diagonal term

D_jj = int( diff(phi_jl,x)^2 ,x ,-hij, 0) +...
       int( diff(phi_jr,x)^2 ,x ,0 ,hjk);               % = 1/hij + 1/hjk

%% overlap 

% O_ij = \int \phi_{i} \phi_{j}

% off diagonal terms

O_ij = int( phi_i*phi_jl, x, -hij, 0);                  % = hij/6

O_jk = int( phi_k*phi_jr ,x ,0 ,hjk);                   % = hjk/6

% diagonal term

O_jj = int( phi_jl^2 ,x ,-hij, 0) +...
       int( phi_jr^2 ,x ,0 ,hjk);                       % = hij/3 + hjk/3


%% potential

% assume the potential is piecewise linear

% offdiagonal terms
V_ij = int( phi_i*phi_jl*((vj-vi)*x/hij+vj),x,-hij,0);  % = hij*(vi+vj)/12

V_jk = int( phi_k*phi_jr*((vk-vj)*x/hjk+vj),x,0,hjk);   % = hjk*(vj+vk)/12

% diagonal term

V_jj = int( phi_jl^2*((vj-vi)*x/hij+vj) ,x ,-hij, 0) +...
       int( phi_jr^2*((vk-vj)*x/hjk+vj) ,x ,0 ,hjk);    % = hij*(vi + 3*vj)/12 + hjk*(3*vj + vk)/12

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Now for something completely different        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
%% Integrals that appear in two problem along boundary

% assume the vibrational wave function is piecewise linear

% B_j = \int \phi_{j}*wf 

B_i = int( phi_jl*((vj-vi)*x/hij+vj),x,-hij,0) +...
      int( phi_jr*((vk-vj)*x/hjk+vj),x,0,hjk);          % = (hij*(vi + 2*vj))/6 + (hjk*(2*vj + vk))/6





