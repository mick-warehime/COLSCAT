% find the reaction path of a given 2D co-linear scattering potential
% using the string method of weinan, weiqing and vanden eijnden
% j chemical phsyics 126 164103 2007

function [phi,vphi,rbar,vbar] = mepath(x,y,v)

%% set parameters
dt = 50;         % time step of steepest descent
tol = 1e-6;      % error ( max displacement from one time step to next)
maxit = 2000;    % second exit condition if number of steps is too high
n = 201;         % number of points in the reaction path
freq = 20;       % number of its between re-spline points to even arc length
plotopt = 1;     % whether or not to plot
vcont = 0.018;    % initial potential contour line
exitmsg = 1;
% conversion from hartree to eV
econv = 27.211396132;

%%

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
xmin = min(min(x)); xmax = max(max(x));
ymin = min(min(y)); ymax = max(max(y));

% initial phi is simply a contour line along the repulsive wall
phil  = contour(x,y,v,[vcont vcont]);
nphil = phil(2,1);
phil  = phil(:,2:nphil);

% make sure there are only n points in the initial phi by splining
pvec = 1:(nphil-2)/(n-1):(nphil-1);
phix  = spline(1:(nphil-1),phil(1,:),pvec);
phiy  = spline(1:(nphil-1),phil(2,:),pvec);
phi = spline_phi([phix',phiy'],n);

% % intial string (ensures we are in both basins of attraction)
% phi = [(xmin:(xmax-xmin)/(n-1):xmax)' (ymax:-(ymax-ymin)/(n-1):ymin)'];

% calculate the gradients
[gx,gy] = gradient(v/econv);

% time counter
count=0;

% set initial values of the potential integral and the error
vint = 0;
err = 10;

% move string 'down hill' until the potential stops changing
while err>tol && count<maxit
    
    % get the potential along the line
    v_phi = interp2(x,y,v,phi(:,1),phi(:,2));
    
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
        vint = trapz(interp2(x,y,v,phi(:,1),phi(:,2)));
        
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
vphi = interp2(x,y,v,phi(:,1),phi(:,2));

% spline to determine barrier location and height
ind = 1:length(vphi);
barin = fminsearch(@(r) -spline(ind,vphi,r),length(ind)/2);

% store relevant path info
rbar = [spline(ind,phi(:,1),barin) spline(ind,phi(:,2),barin)];
vbar = spline(ind,vphi,barin);


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
