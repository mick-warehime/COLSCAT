function d=dimp_hyp(p,par,nit,alpha)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

if nargin<3, nit=20; end
if nargin<4, alpha=0.1; end

x0=p(:,1);
y0=p(:,2);
x=x0;
y=y0;
for it=1:nit
  [cf cfx cfxx cfy cfyy] = hypdif(p,par);

  F1=cf;
  F2=(x-x0).*cfy-(y-y0).*cfx;
  J11=cfx;
  J12=cfy;
  J21=cfy-(y-y0).*cfxx;
  J22=-cfx+(x-x0).*cfyy;
  
  detJ=J11.*J22-J12.*J21;
  detJ(detJ==0)=inf;
  
  x=x-alpha*(J22.*F1-J21.*F2)./detJ;
  y=y-alpha*(-J12.*F1+J11.*F2)./detJ;
end

d=sqrt((x-x0).^2+(y-y0).^2).*sign(hypdif([x0,y0],par));
