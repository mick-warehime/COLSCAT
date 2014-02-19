% return the divergence of a function defined on a triangulated grid

function div = tridivergence(t,x,y,z)

trisize = size(t,2);

if trisize==3    
    [zxx,junk] = trigradientP1(x,y,z(:,1),t);
    [junk,zyy] = trigradientP1(x,y,z(:,2),t);
elseif trisize==6
    [zxx,junk] = trigradientP2(x,y,z(:,1),t);
    [junk,zyy] = trigradientP2(x,y,z(:,2),t);
end
div = zxx + zyy;

