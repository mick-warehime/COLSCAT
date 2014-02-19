% the generalized morse potential

function vgm = genmorse(r,par)

% par = [a re de gamma]
% v_generalized morse = de*{ (gamma/a)*[exp(-2a*(r-re))-exp(-a(r-re))]+1-exp(-gamma(r-re) }

vgm = par(3)*( (par(4)/par(1)) *(exp(-2*par(1)*(r-par(2)))-exp(-par(1)*(r-par(2))))+1-exp(-par(4)*(r-par(2))));

% plot(r,vgm,'ko')
% drawnow