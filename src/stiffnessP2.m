
function A = stiffnessP2(p,t)

% points in triangulation
x = p(:,1); y = p(:,2);
x = x(t); y = y(t);

% transformation from msj-coordinates to standard triangle coordinates
a = y(:,3)-y(:,1); b = y(:,1)-y(:,2);
c = x(:,1)-x(:,3); d = x(:,2)-x(:,1); 

abcd = [a.^2 a.*b b.^2 c.^2 c.*d d.^2];


% coef comes form the stiff matrix elements in p2integrals.m
coef = [3     6     3     3     6     3;
        1     1     0     1     1     0;
        0     1     1     0     1     1; 
       -4    -4     0    -4    -4     0
        0     0     0     0     0     0;
        0    -4    -4     0    -4    -4;
        1     1     0     1     1     0;
        3     0     0     3     0     0;
        0    -1     0     0    -1     0;
       -4    -4     0    -4    -4     0;
        0     4     0     0     4     0;
        0     0     0     0     0     0;
        0     1     1     0     1     1;
        0    -1     0     0    -1     0;
        0     0     3     0     0     3;
        0     0     0     0     0     0;
        0     4     0     0     4     0;
        0    -4    -4     0    -4    -4;
       -4    -4     0    -4    -4     0;
       -4    -4     0    -4    -4     0;
        0     0     0     0     0     0;
        8     8     8     8     8     8;
        0    -8    -8     0    -8    -8;
        0     8     0     0     8     0;
        0     0     0     0     0     0;
        0     4     0     0     4     0;
        0     4     0     0     4     0;
        0    -8    -8     0    -8    -8;
        8     8     8     8     8     8;
       -8    -8     0    -8    -8     0;
        0    -4    -4     0    -4    -4;
        0     0     0     0     0     0;
        0    -4    -4     0    -4    -4;
        0     8     0     0     8     0;
       -8    -8     0    -8    -8     0;
        8     8     8     8     8     8]/6;

 A = zeros(length(x),36);
 for j = 1:36     
     A(:,j) = sum(bsxfun(@times,abcd,coef(j,:)),2);
 end
