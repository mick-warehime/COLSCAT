% convert to bond coordinates

function U = tobond(R,m)

M = mass(m);

% convert from (Ra,ra) to (U1,U2)
U(:,2) = R(:,2)*M.lama;
U(:,1) = R(:,1)/M.lama-U(:,2)*m(3)/(m(2)+m(3));
end