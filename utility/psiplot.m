% plot a multistate wave function

function psiplot(mesh,psi)

t = mesh.t;
p = mesh.p;
pb = [mesh.pb; mesh.pb(1,:)];

figure(1)
set(gcf,'position',[447,526,1140,483])


subplot(1,2,1)
title('\psi_1')
tricontour(p,t,abs(psi(:,1)).^2,0.05:0.05:1)
hold on
plot(pb(:,1),pb(:,2),'k')
xlabel('R_a')
ylabel('r_a')

subplot(1,2,2)
title('\psi_2')
tricontour(p,t,abs(psi(:,2)).^2,0.05:0.05:1)
hold on
plot(pb(:,1),pb(:,2),'k')
xlabel('R_a')
ylabel('r_a')

% 
% subplot(1,2,1)
% title('\psi_1')
% trisurf(t,p(:,1),p(:,2),abs(psi(:,1)).^2,'edgecolor','none','facecolor','interp')
% view(35,50)
% 
% subplot(1,2,2)
% title('\psi_2')
% trisurf(t,p(:,1),p(:,2),abs(psi(:,2)).^2,'edgecolor','none','facecolor','interp')
% view(35,50)