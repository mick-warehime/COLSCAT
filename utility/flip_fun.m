% tanh flipping function

function s = flip_fun(x, x_0, alpha, beta)

% returns the value s = alpha/2*(1-tanh(beta*(theta-theta_0)))

s = alpha/2*(1-tanh(beta*(x-x_0)/pi));