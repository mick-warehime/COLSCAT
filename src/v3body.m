function v3 = v3body(x,y,mab,mbc,D)
morseab = genmorse(x,mab);
morsebc = genmorse(y,mbc);
v3 = zeros(size(x,1),size(D,2));
for j=1:size(D,2);
    v3(:,j) = morseab+morsebc+threebody(D(:,j)',[x y]);
end

% force the potential to be asymptotically vdif == 0
sx = flip_fun(x,4,1,6.5);
sy = flip_fun(y,4,1,6.5);
s = sx.*sy;
v3(:,1) = v3(:,1).*s+(1-s).*v3(:,2); 
