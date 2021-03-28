function [xmoy,pmoy,incert_x,incert_p,E,t,Pgauche,Pdroite,V,x,psi2]=AnalyseSchroedinger(fichier,omega)

%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fichier = 'output';
data = load([fichier,'_obs.dat']);
t = data(:,1);
Pgauche = data(:,2);
Pdroite = data(:,3);
E = data(:,4);
xmoy = data(:,5);
x2moy = data(:,6);
pmoy = data(:,7);
p2moy = data(:,8);
incert_x=data(:,9);
incert_p=data(:,10);
data = load([fichier,'_pot.dat']);
x = data(:,1);
V = data(:,2);
psi2 = load([fichier,'_psi2.dat']);

%% Figures %%
%%%%%%%%%%%%%
% figure('Name',['Analyse de ' fichier])
% subplot(2,2,1)
% plot(x,V), hold on
% plot(x([1,end]),E(1)*ones(1,2),'--')
% grid on
% xlabel('x')
% ylabel('V')
% legend('V(x)','E','Location','Best')
% 
% subplot(2,2,2)
% [X,T] = meshgrid(x,t);
% pcolor(X,T,psi2)
% shading interp
% colormap jet
% c = colorbar;
% xlabel('x [m]')
% ylabel('t [s]')
% ylabel(c,'|\psi|^2')
% 
% subplot(2,2,3)
% plot(t,Pgauche,t,Pdroite,t,Pgauche+Pdroite)
% grid on
% xlabel('t')
% ylabel('P')
% legend('P_{x<0}','P_{x>0}','P_{tot}','Location','Best')
% 
% subplot(2,2,4)
% plot(t,E)
% grid on
% xlabel('t')
% ylabel('E')
% 
% X1=[0,1000];
% Y1=[0.5,0.5];
% 
% figure (100)
% hold on
% plot(t,incert_x.*incert_p,'.')
% plot(X1,Y1,'k--')
% set(gca,'FontSize',16)
% xlabel('t')
% ylabel('<\Delta x>(t)\cdot<\Delta p>(t)')
% 
% X=[0,1000];
% Y=[sqrt(1/(2*omega(end))),sqrt((1/(2*omega(end))))];
% [peaks, loc]=findpeaks(incert_x);
% t_peaks=t(loc);
% 
% 
% figure (11)
% hold on
% plot(t,incert_x)
% plot(X,Y,'--')
% plot(t_peaks,peaks,'g')
% set(gca,'FontSize',16)
% xlabel('t')
% ylabel('<\Deltax>')
% legend({'$\left<\Delta x\right>$','$\sqrt{\frac{\hbar}{2m\omega}}$','$f(x)=7\cdot10^{-4}x+5$'},'Interpreter','latex')

save('data')