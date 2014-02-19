function v=potfh2_col(rfh,rhh)
%     Stark-Werner F+H2 potential energy surface

%     r(1) = H-H distance in bohr
%     r(2) = F-H distance in bohr
%     r(3) = H-F distance in bohr
 

% v = potential in hartree from bottom of asymptotic F+H2 valley

% initialize fh2 parameters
[a,b,c,d,ma,p] = datfh2_col();

% determine bond coordinates for collinear geometry
rb = [rhh(:)  rfh(:) rfh(:)+rhh(:)];

%     (Modified) Aguado-Paniagua type functions:

x = min(rb(:,2:3),[],2);
y = rb(:,1);
z = max(rb(:,2:3),[],2);
b1 = a(ma-5);
b2 = a(ma-4);
b3 = a(ma-3);
x0 = a(ma-2);
y0 = a(ma-1);
z0 = a(ma);

B1 = b1*bsxfun(@times,b(1:ma-6)',x-x0);
B2 = b2*bsxfun(@times,c(1:ma-6)',y-y0);
B3 = b3*bsxfun(@times,d(1:ma-6)',z-z0);
fex=exp(-(B1+B2+B3));

X = bsxfun(@power,x,b(1:ma-6)');
Y = bsxfun(@power,y,c(1:ma-6)');
Z = bsxfun(@power,z,d(1:ma-6)');
fxy=X.*Y.*Z.*fex;

fit=sum(bsxfun(@times,a(1:ma-6)',fxy),2);


%      xr=x-p(3);
xr=[x z]-p(3);
yr=y-p(9);
%     zr=z-p(3);
fx=exp(-p(2)*xr);
fy=exp(-p(8)*yr);
%      fz=exp(-p(2)*zr);
xv=-p(1)*polyval([p(5) p(4) p(2) 1],xr).*fx+p(6);
yval=-p(7)*polyval([p(11) p(10) p(8) 1],yr).*fy+p(12);

%       xval=-p(1)*(1+p(2)*xr+p(4)*xr^2+p(5)*xr**3)*fx+p(6);
%       yval=-p(7)*(1+p(8)*yr+p(10)*yr^2+p(11)*yr**3)*fy+p(12);
%       zval=-p(1)*(1+p(2)*zr+p(4)*zr**2+p(5)*zr**3)*fz+p(6);
%      vagpan=fit+xval+yval+zval;

% add the h-h potential
v=fit+sum(xv,2)+yval;


return
