% plot the mesh boundary
function b = meshboundary(mesh)


pr_fi = protate(mesh.pf,mesh.pi(5));
pr_fo = protate(mesh.pf,mesh.po(5));

% par = [a b xc yc]
n = 0:99;

xnot_i = mesh.pi(3)+mesh.pi(1);
xi_top = fliplr(xnot_i+n*(pr_fi(2,1)-xnot_i)/n(end));
xi_bot = xnot_i+n*(pr_fi(3,1)-xnot_i)/n(end);
yi_top = mesh.pi(2)*sqrt( ((xi_top-mesh.pi(3))/mesh.pi(1)).^2-1)+mesh.pi(4);
yi_bot = -mesh.pi(2)*sqrt( ((xi_bot-mesh.pi(3))/mesh.pi(1)).^2-1)+mesh.pi(4);


xnot_o = mesh.po(3)+mesh.po(1);
xo_top = fliplr(xnot_o+n*(pr_fo(1,1)-xnot_o)/n(end));
xo_bot = xnot_o+n*(pr_fo(4,1)-xnot_o)/n(end);
yo_top = mesh.po(2)*sqrt( ((xo_top-mesh.po(3))/mesh.po(1)).^2-1)+mesh.po(4);
yo_bot = -mesh.po(2)*sqrt( ((xo_bot-mesh.po(3))/mesh.po(1)).^2-1)+mesh.po(4);

pi = protate([[xi_top,xi_bot]',[yi_top,yi_bot]'],-mesh.pi(5));
po = protate([[xo_top,xo_bot]',[yo_top,yo_bot]'],-mesh.po(5));

b = [pi; mesh.pf(3:4,:); flipud(po); mesh.pf(1:2,:)];


