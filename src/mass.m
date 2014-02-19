% returns the mass structure with relevant masses calculated
function M = mass(m)

% atomic masses
M.ma = m(1);
M.mb = m(2);
M.mc = m(3);

% total mass
M.M = sum(m);

% skew angle
M.alpha = alph(m);

% reduced masses
M.mu   = sqrt(m(1)*m(2)*m(3)/M.M);       % system
M.muab = m(1)*m(2)/(m(1)+m(2));          % product  diatomic
M.mubc = m(2)*m(3)/(m(2)+m(3));          % reactant diatomic
M.muc  = m(3)*(m(1)+m(2))/M.M;           % product  arrangement
M.mua  = m(1)*(m(2)+m(3))/M.M;           % reactant arrangement

% lambda factors
M.lama = sqrt(M.mu/M.mubc);              % reactant
M.lamc = sqrt(M.mu/M.muab);              % product
 
end

