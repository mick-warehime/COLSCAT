function [w,theta] = normalmodes(m,rbar,h,vinput)

% mass matrix of KE operator
M = [m(1)*(m(2)+m(3)) m(1)*m(3); m(1)*m(3) m(3)*(m(1)+m(2))]/sum(m);

% diagonalize M
[C,L] = eig(M);

% build D matrix (divide columns by the square root of the eigenvalues)
D = C*diag((1./sqrt(diag(L))));

% determine the potential in small grid around the barrier
R = bsxfun(@plus,rbar,h*[1 1; 1 0; 1 -1; 0 1; 0 0; 0 -1; -1 1; -1 0; -1 -1]);
V = feval(vinput,R(:,1),R(:,2),1,[]);

% calculate the hessian at the barrier
hxx = (sum(V([1:3:9 3:3:9]))-2*sum(V(2:3:8)))/3;
hyy = (sum(V(1:3))-2*sum(V(4:6))+sum(V(7:9)))/3;
hxy = (sum(V([1 9]))-sum(V([3 7])))/4;
H = [hxx hxy; hxy hyy]/h^2;

% rotate the potential matrix. 
F = D'*H*D;

% calculate the vibrational frequencies
[K,k] = eig(F);

% normal mode frequence
w = sqrt(max(k(:)));

% determine the rotation to the normal coordinates
vec=[ 0 0 -1 1;-1 1 0 0]*.5;
theta = K*C*vec;


