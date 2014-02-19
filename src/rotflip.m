% rotate the system by alpha degrees and then invert the y axis

function [x y] = rotflip(x,y,alpha)

% R = [cos(alpha) sin(alpha); sin(alpha) -cos(alpha);
    
x = cos(alpha)*x + sin(alpha)*y;
y = sin(alpha)*x - sin(alpha)*y;
