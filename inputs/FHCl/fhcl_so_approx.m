% an approximation to the fhcl so parameter

function so = fhcl_so_approx(x,y,b)


% points on the dividing line (picked using ginput)
% xd = [2.0025 1.9983; 2.6062 3.0044; 3.3027 3.9950];
xd = [1.9438    2.0852
    2.3245    2.4213
    2.6906    3.1055
    3.0127    3.5616
    3.3421    4.0898
    3.9205    4.7620];
coef  = polyfit(xd(:,1),xd(:,2),1);

% calculate the signed distance between each point and the line
d = (-y+coef(1)*x+coef(2))/sqrt(coef(1)^2+1);

% create the flipping function as a function of the signed distance
fl = flip_fun(d, 0, 1, b);

% approx asytmptotic value of the so value for f in hartree
sof = 129.5/219474.63;

% approx asytmptotic value of the so value for cl in hartree
socl = 273/219474.63; 

% create the desired flipped output
so = fl*socl+(1-fl)*sof;

% % approximate vdiff
% xc = [2.35 2.55];
% vdif = .1*exp((-(x-xc(1)).^2-(y-xc(2)).^2)*2);