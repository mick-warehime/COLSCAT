function N = hypfit(p,par)

pr = protate(p,par(5));
 
% model the x values using a hyperbolic fit ((x-xc)/a)^2-(y-yc)^2/b^2=1
% par = [a b xc yc]
xfit = par(1)*sqrt(1+(pr(:,2)-par(4)).^2/par(2)^2)+par(3);

% return the normal (i.e. set ans = 2 norm of zfit-z)
N = norm(abs(xfit-pr(:,1))./pr(:,1).^2 );


% %% plotting diagnostics
% plot(pr(:,1),pr(:,2),'o-',xfit,pr(:,2),'ro-')
% drawnow


