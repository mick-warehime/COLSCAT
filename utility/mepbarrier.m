% find the reaction path of a given 2D co-linear scattering potential
% using the string method of weinan, weiqing and vanden eijnden
% j chemical phsyics 126 164103 2007

function [rbar,vbar] = mepbarrier(rbar0,theta,dr,vinput)

%% set parameters
nsearch = 100;    % number of points along each search direction
dt = 100;         % time step of steepest descent
tol = 1e-6;      % error ( max displacement from one time step to next)
maxit = 100;    % second exit condition if number of steps is too high
n = 101;         % number of points in the reaction path
freq = 20;       % number of its between re-spline points to even arc length
plotopt = 1;     % whether or not to plot
exitmsg = 1;     


%%
eps = 0.8;
xvec = rbar0(1)+dr*cos(theta*pi/180)*[-1 1]*eps;
yvec = rbar0(2)+dr*sin(theta*pi/180)*[1 -1]*eps;

% create search meshgrid
[x,y] = meshgrid((rbar0(1)-dr):dr/(nsearch-1):(rbar0(1)+dr),(rbar0(2)-dr):dr/(nsearch-1):(rbar0(2)+dr));
v = reshape(feval(vinput,x,y),size(x));
v = v-interp2(x,y,v,rbar0(1),rbar0(2));

%intialize search string in the theta direction
phi =  [linspace(xvec(1),xvec(2),n); linspace(yvec(1),yvec(2),n)]';

% plot contour if we are plotting
if plotopt
    figure
    hold on
    
    % potential contours values
    dv = abs(max(v(:)) - min(v(:)))/500;
    vcons = min(v(:)):dv:30*dv;
    
    % plot contour
    contour(x,y,v,vcons)
    
end

% range of mesh grid
xmin = rbar0(1)-dr; xmax = rbar0(1)+dr;
ymin = rbar0(2)-dr; ymax = rbar0(2)+dr;

% calculate the gradients
[gx,gy] = gradient(v);

% time counter
count=0;

% set initial values of the potential integral and the error
vint = 0;
err = 10;

% move string 'down hill' until the potential stops changing
while err>tol && count<maxit
  
    % get the potential along the line
    v_phi = feval(vinput,phi(:,1),phi(:,2));    
    
    % calculate normal vector ( the string only moves in the normal dir)
    normal = getnorm(phi,v_phi);
    
    % calculate force vector fx = -dv/dx, fy = -dv/dy, along the string
    fx = interp2(x,y,gx,phi(:,1),phi(:,2));
    fy = interp2(x,y,gy,phi(:,1),phi(:,2));
    
    % calculate the displacement
    dphi = [fx.*normal(:,1) fy.*normal(:,2) ];
    
    % move the points
    phi = phi-dt*dphi;

    % move points back into boundary in case they move outside
    phi(phi(:,1)>xmax,1) = xmax; phi(phi(:,1)<xmin,1) = xmin;
    phi(phi(:,2)>ymax,2) = ymax; phi(phi(:,2)<ymin,2) = ymin;
    
    % interpolate to get evenly spaced points along the line & compute error
    if mod(count,freq)==0
        
        % respline the points to be evenly spaced along the path
        phi = spline_phi(phi,n);
        
        % store old potential integral and new potential integral
        vint_old = vint;
        vint = trapz(feval(vinput,phi(:,1),phi(:,2)));
        
        % error depends on the change in the integral of the potential
        err = abs(vint-vint_old);
        
    end
    
    count = count + 1;
    
    % plot every ten frames
    if count == 1 && plotopt
        ht = plot(phi(:,1),phi(:,2),'ko');
        title(['time step: ' num2str(count), '. error: ', num2str(err, '%4.2e \n')])
        drawnow
    elseif mod(count,10)==0  && plotopt
        set(ht,'XData', phi(:,1));
        set(ht,'YData', phi(:,2));
        drawnow
        title(['time step: ' num2str(count), '. error: ', num2str(err, '%4.2e \n')])
    end
 
 end
close


% potential energy along reaction path
vphi = feval(vinput,phi(:,1),phi(:,2));

% spline to determine barrier location and height
ind = 1:length(vphi);
barin = fminsearch(@(r) -spline(ind,vphi,r),length(ind)/2);

% store relevant path info
rbar = [spline(ind,phi(:,1),barin) spline(ind,phi(:,2),barin)];
vbar = feval(vinput,rbar(1),rbar(2));

% if an exit message is wanted, give them what they want
if exitmsg
    % quick text to describe exit strategy
    if err<tol
        disp('Error less than tolerance. Exit.')
    elseif count>maxit
        disp('Reached max iterations. Abort.')
    end
end

end

% compute the direction of the normal to the string at each point
% using finite differences based on the rules on page 164103-8 of ref
function normal = getnorm(phi,v_phi)
vin = v_phi(2:end)>v_phi(1:end-1);
vin = [0; vin];

lv = length(vin);

tangent = zeros(lv,2);

for j=1:length(vin)
    
    
    if j==1
        tangent(j,:) = (phi(2,:)-phi(1,:))/norm(phi(2,:)-phi(1,:));
    elseif j==length(vin)
        tangent(j,:) = (phi(end,:)-phi((end-1),:))/norm(phi(end,:)-phi((end-1),:));
    else
        % determine the proper direction to calculate the tangent vector
        if vin(j) && vin(j-1)
            tangent(j,:) = (phi(j+1,:)-phi(j,:))/norm(phi(j+1,:)-phi(j,:));
        elseif ~vin(j) && ~vin(j-1)
            tangent(j,:) = (phi(j,:)-phi(j-1,:))/norm(phi(j,:)-phi(j-1,:));
        elseif (~vin(j)  && vin(j-1)) || (vin(j) && ~vin(j-1))
            tangent(j,:) = (phi(j+1,:)-phi(j-1,:))/norm(phi(j+1,:)-phi(j-1,:));
        end
    end
    
    
end
% return the vector normal to the tangent
normal = [-tangent(:,2) tangent(:,1)];
end

% spline the points to ensure they are equally spaced along the string
function phi = spline_phi(phi,n)

% finer grid in the x coordinate
xx = 1:(n-1)/1000:n;

% spline to get the arc lengths
sx = spline(1:n,phi(:,1),xx);
sy = spline(1:n,phi(:,2),xx);

dsx = sx(2:end)-sx(1:(end-1));
dsy = sy(2:end)-sy(1:(end-1));

% create an arclength matrix
svec = sqrt(dsx.^2+dsy.^2);
smat = repmat(svec',1,length(dsx));

% approximation of the arc length
s = trapz(triu(smat));

% vector with n points and equal arc length between each point
ds = s(end)/(n-1);
sn = 0:ds:s(end);

% spline to get new indices
in = spline([0 s],xx,sn);

% now get x and y values
phi(:,1) = spline(1:n,phi(:,1),in);
phi(:,2) = spline(1:n,phi(:,2),in);


end
