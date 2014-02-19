
void

%% 1D P2 integrals
x = sym('x','real');
dx = sym('dx','real');
syms v1 v2 v3

C2_1D = [0 0 1; dx^2 dx 1; 4*dx^2 2*dx 1]\eye(3);
C2_v1D = [0 0 1; dx^2 dx 1; 4*dx^2 2*dx 1]\[v1;v2;v3];

phi2_1d = C2_1D.'*[x^2;x;1];


for i=1:3
    for j=1:3
        O_1D2(i,j) = int(phi2_1d(i)*phi2_1d(j),x,0,2*dx);
        T_1D2(i,j) = int(diff(phi2_1d(i),x)*diff(phi2_1d(j),x),x,0,2*dx);
        for l = 1:3
            V_1D2(i,j,l) = int(phi2_1d(i)*phi2_1d(l)*phi2_1d(j),x,0,2*dx);
        end
    end
end

%% 1D P1 Integrals

x = sym('x','real');
dx = sym('dx','real');
syms v1 v2

C1_1D  = [0 1; dx 1]\eye(2);

phi1_1d = C1_1D.'*[x;1];

for i=1:2
    for j=1:2
        
        O_1D1(i,j) = int(phi1_1d(i)*phi1_1d(j),x,0,dx);
        T_1D1(i,j) = int(diff(phi1_1d(i),x)*diff(phi1_1d(j),x),x,0,dx);
        for l = 1:2
            V_1D1(i,j,l) = int(phi1_1d(i)*phi1_1d(j)*phi1_1d(l),x,0,dx);
        end
    end
end



