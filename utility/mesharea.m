% determine the approximate area of the mesh

function area = mesharea(mesh,n)
dbox = diff(mesh.hbox);

x = rand(n,1)*dbox(1)+mesh.hbox(1,1);
y = rand(n,1)*dbox(1)+mesh.hbox(1,2);

d = dmesh([x(:),y(:)],mesh);
in = d<0;

area = (length(find(in))/length(x(:)))*prod(dbox);
