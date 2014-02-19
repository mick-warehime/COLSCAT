% vector valued function of the three body fit
function v3 = threebody(par,R)

M = (sqrt(9-8*(1-(length(par)-2)))-3)/2;

beta = par(1:2);
D = par(3:end);

r1 = R(:,1).*exp(-R(:,1)*beta(1));
r2 = R(:,2).*exp(-R(:,2)*beta(2));

[a,b] = meshgrid(0:M,0:M);
c = [b(:) a(:)];
d = c(sum(c,2)<=M,:);

Phi = bsxfun(@power,r1,d(:,1)').*bsxfun(@power,r2,d(:,2)');

v3 = sum(bsxfun(@times,D,Phi),2);