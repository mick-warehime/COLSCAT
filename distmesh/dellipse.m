function d=dellipse(p,pb,a,b,alpha,xc,ps)

%   Copyright (C) 2012 MICK  WAREHIME

% set the intersectino points to the origin
ps1 = pshift(p,-ps(1,1),-ps(1,2));
ps2 = pshift(p,-ps(2,1),-ps(2,2));

% rotate by alpha/2
pr1 = protate(ps1,alpha/2);
pr2 = protate(ps2,alpha/2);

% calculate distance to ellipse
de1 = sqrt( ((pr1(:,1)-xc(1))/a(1)).^2+(pr1(:,2)/b(1)).^2)-1;
de2 = sqrt( ((pr2(:,1)-xc(2))/a(2)).^2+(pr2(:,2)/b(2)).^2)-1;

% calculate distance to polygon
dp = dpoly(p,pb);

% we want the difference between the  union of de2 and dp and de1
d = ddiff(dunion(de2,dp),de1);



