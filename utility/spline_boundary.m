function phi = spline_boundary(phi,n,eps)
   
% tolerance for percent difference in slope
mtol = 0.025;

% phi = even_spline(phi,10*n);
% ------------------------------------------------------------
% remove points that fit a straight line with in mtol

% determine the slope of the line between pts 1 and 2
coef_1 = polyfit(phi(1,1:2),phi(2,1:2),1);

for j=3:length(phi)
   %
    % calculate the slope betwweens pts 1 and j
    coef_j = polyfit([phi(1,1) phi(1,j)],[phi(2,1) phi(2,j)],1);
    
    % if the slope between pts 1 and j is greater than mtol break
    if abs((coef_1(2)-coef_j(2))/coef_1(2))>mtol       
        in1 = j-1;
        break
    end
    
end

% determine the slope of the line between pts end and end-1
coef_end = polyfit(phi(1,(end-1):end),phi(2,(end-1):end),1);

for j=2:length(phi)
   
    % calculate the slope betwweens pts end and (end-j)
    coef_j = polyfit([phi(1,end) phi(1,(end-j))],[phi(2,end) phi(2,(end-j))],1);    
    
    % if the slope between pts end and (end-j) is greater than mtol break
    if abs((coef_end(2)-coef_j(2))/coef_end(2))>mtol       
        inend = length(phi)-j+1;
        break
    end
    
end

% remove points that are along a straight line


phi(:,[2:in1 inend:(end-1)]) = [];

% % spline to put more points around the curved part of the boundary
% lp = length(phi);
% [~, in] = min(phi(1,:));
% in = in+eps;
% d = [0 erf((1:(n-1))/(n/2)) 1 ];
% sin = [1 2+(in-2)*d in*( (1-d((end-1):-1:2))+1) lp];
% 
% phi = spline(1:length(phi),phi,sin);

end


function phiN = even_spline(phi,n) 

np = length(phi);

% finer grid in the x coordinate
xx = 1:(np-1)/1000:np;

% spline to get the arc lengths
sx = spline(1:np,phi(1,:),xx);
sy = spline(1:np,phi(2,:),xx);

dsx = sx(2:end)-sx(1:(end-1));
dsy = sy(2:end)-sy(1:(end-1));

% create an arclength matrix
svec = sqrt(dsx.^2+dsy.^2);
smat = repmat(svec',1,length(dsx));

% approximation of the arc length
s = trapz(triu(smat));

% vector with n points and equal arc length between each point
ds = max(s)/(n-1);
sn = 0:ds:max(s);

% spline to get new indices
in = spline([0 s],xx,sn);

% now get x and y values
phiN(1,:) = spline(1:np,phi(1,:),in);
phiN(2,:) = spline(1:np,phi(2,:),in);

end