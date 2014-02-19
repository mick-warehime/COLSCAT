% return the intersection points between the hyperbolic curves and the
% boundaries

function pfix = hypline(pi,po,rmax,pmax,alph)

% lines for product/reactant boundaries in (Ra,ra) coordinates

% reactant boundary x = rmax in (Ra,ra)
pfix(4,:) = protate(hlint(tan(po(5)),rmax/cos(po(5)),po),-po(5));
pfix(3,:) = protate(hlint(tan(pi(5)),rmax/cos(pi(5)),pi),-pi(5)); 


% product boundary y = -1/tan(alph)x + pmax/sin(alph) in (Ra,ra)
ai = cos(pi(5))*sin(alph)-sin(pi(5))*cos(alph); bi = cos(pi(5))*cos(alph)+sin(pi(5))*sin(alph);
ao = cos(po(5))*sin(alph)-sin(po(5))*cos(alph); bo = cos(po(5))*cos(alph)+sin(po(5))*sin(alph);

pfix(2,:) = protate(hlint(-ai/bi,pmax/bi,pi),-pi(5));
pfix(1,:) = protate(hlint(-ao/bo,pmax/bo,po),-po(5));
end

function x = hlint(M,B,p)

% find the intersection between the line x = My+B and the hyperbola

% ((x-xc)/a)^2-((y-yc)/b)^2 -1 =0;
% ((x-p(3))/p(1)).^2-((y-p(4))/p(2)).^2-1 =0;
%x = p(1)*sqrt(((y-p(4))/p(2)).^2+1)+p(3);
a=(M^2/p(1)^2 - 1/p(2)^2);
b=((2*p(4))/p(2)^2 - (2*M*p(3))/p(1)^2 + (2*B*M)/p(1)^2); 
c=B^2/p(1)^2 + p(3)^2/p(1)^2 - p(4)^2/p(2)^2 - (2*B*p(3))/p(1)^2 -1;

xp(2) = (-b+sqrt(b^2-4*a*c))/(2*a);
xm(2) = (-b-sqrt(b^2-4*a*c))/(2*a);
xp(1)  = M*xp(2)+B;
xm(1)  = M*xm(2)+B;

if M<0
    x = xm;
else
    x = xp;
end

end

