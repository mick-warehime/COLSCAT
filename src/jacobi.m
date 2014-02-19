% determine the reactant and product jacobi coordinates
function jac = jacobi(p,bnd,M)
jac.Ra = p(bnd.ba,1)/M.lama;
jac.ra = p(bnd.ba,2)*M.lama;
jac.Rc = (cos(M.alpha)*p(bnd.bc,1)+sin(M.alpha)*p(bnd.bc,2))/M.lamc;
jac.rc = (sin(M.alpha)*p(bnd.bc,1)-cos(M.alpha)*p(bnd.bc,2))*M.lamc;