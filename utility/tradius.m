function k = tradius(x,y,t)

x = x(t);   y = y(t);

a = sqrt( (x(:,1)-x(:,2)).^2 + (y(:,1)-y(:,2)).^2);
b = sqrt( (x(:,2)-x(:,3)).^2 + (y(:,2)-y(:,3)).^2);
c = sqrt( (x(:,1)-x(:,3)).^2 + (y(:,1)-y(:,3)).^2);

s = (a+b+c)/2;
area = sqrt(s.*(s-a).*(s-b).*(s-c));

k = a.*b.*c./(4*area);

 