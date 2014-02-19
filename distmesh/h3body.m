function h = h3body(p,pfix,pi,po,vfun,m)

pb = tobond(p,m);

v = vfun(pb(:,1),pb(:,2));

h = 2-exp(-(v(:,1)+.025)).^8;