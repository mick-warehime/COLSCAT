% fit an ellipse to the corner of each boundary

function N = ellfit(x,y,par)
% par = [a b xc yc theta]

% rotate by theta
p = protate([x; y]',par(5));

% shift by xc yc
p = pshift(p,-par(3),-par(4));

% scale by a and b
p = [p(:,1)/par(1) p(:,2)/par(2)];

% model the x values using a elliptical fit ((x-xc)/a)^2+((y-yc)/b)^2=1
xfit = -sqrt(1-p(:,2).^2);

% return the normal (i.e. set ans = 2 norm of zfit-z)
N = sum( abs( p(:,1) -xfit).*abs(xr));

plot(p(:,1),p(:,2),'ro',xfit,p(:,2),'ko')
drawnow
