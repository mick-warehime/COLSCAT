% calculate the curl of a vector field defined above a triangulated mesh

function curlz = tricurl(t,x,y,fx,fy)

% determine if the triangulation is linear or quadratic
[~,trisize] = size(t);

if trisize==3    
    [junk,fxy] = trigradientP1(x,y,fx,t);
    [fyx,junk] = trigradientP1(x,y,fy,t);
elseif trisize==6
    [junk,fxy] = trigradientP2(x,y,fx,t);
    [fyx,junk] = trigradientP2(x,y,fy,t);
end

curlz = fyx-fxy;